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
from sage.manifolds.chart import Chart
from sage.manifolds.differentiable.real_line import OpenInterval
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
      :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet`
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
    - ``parameters`` -- list of the symbolic expressions used in
      ``equations_rhs`` and ``initial_tangent_vector`` other than the
      coordinates, the velocities and the curve parameter
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

        sage: M = Manifold(3, 'M')
        sage: X.<x1,x2,x3> = M.chart()
        sage: var('t B_0 m q L T')
        (t, B_0, m, q, L, T)
        sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
        sage: D = X.symbolic_velocities(); D
        [Dx1, Dx2, Dx3]
        sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]

    Set the initial conditions::

        sage: p = M.point((0,0,0), name='p')
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((1,0,1))

    Declare an integrated curve and display information relative to it::

        sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
        ....:                        parameters=[B_0, m, q, L, T]); c
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
        d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
        d(x2)/dt = Dx2
        d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
        d(x3)/dt = Dx3
        d(Dx3)/dt = 0
        <BLANKLINE>

    Generate a solution of the system and an interpolation of this
    solution::

        sage: sol = c.solve(step=0.2,
        ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:1},
        ....:         solution_key='carac time 1')
        Performing 4th order Runge-Kutta integration by default...
        Numerical integration completed. Resulting list of points was
         associated with the key 'carac time 1' (if this key already
         referred to a former numerical solution, such a solution was
         erased).
        sage: interp = c.interpolate(solution_key='carac time 1',
        ....:                        interpolation_key='interp 1')
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

        sage: c_plot_2d_1 = c.plot(ambient_coords=[x1, x2],
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
        eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
                               parameters=[B_0, m, q, L, T])
        c.solve(step=0.2,
                parameters_values={B_0:1, m:1, q:1, L:10, T:1},
                solution_key='carac time 1')
        c.interpolate(solution_key='carac time 1',
                      interpolation_key='interp 1')
        c_plot_2d_1 = c.plot(ambient_coords=[x1, x2],
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
        sage: c_plot_3d_100 = c.plot(interpolation_key='interp 100',
        ....:                   thickness=2.5, display_tangent=True,
        ....:                   plot_points=200, plot_points_tangent=10,
        ....:                   scale=0.5, color='green',
        ....:                   color_tangent='orange', verbose=False)
        sage: c_plot_3d_1 = c.plot(interpolation_key='interp 1',
        ....:                   thickness=2.5, display_tangent=True,
        ....:                   plot_points=200, plot_points_tangent=10,
        ....:                   scale=0.5, color='blue',
        ....:                   color_tangent='red', verbose=False)
        sage: viewer3D = 'threejs'
        sage: graph = c_plot_3d_1 + c_plot_3d_100
        sage: graph.show(viewer = viewer3D)

    .. PLOT::

        M = Manifold(3, 'M')
        X = M.chart('x1 x2 x3')
        var('x1 x2 x3 t B_0 m q L T')
        B = B_0*t/T*exp(-(x1**2 + x2**2)/L**2)
        D = X.symbolic_velocities()
        eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
                               parameters=[B_0, m, q, L, T])
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
        c_plot_3d_1 = c.plot(interpolation_key='interp 1',
                        thickness=2.5, display_tangent=True,
                        plot_points=200, plot_points_tangent=10,
                        scale=0.5, color='blue', color_tangent='red',
                        verbose=False)
        c_plot_3d_100 = c.plot(interpolation_key='interp 100',
                            thickness=2.5, display_tangent=True,
                            plot_points=200, plot_points_tangent=10,
                            scale=0.5, color='green',
                            color_tangent='orange', verbose=False)
        graph = c_plot_3d_1 + c_plot_3d_100
        sphinx_plot(graph)

    """

    def __init__(self, parent, equations_rhs, velocities,
                 curve_parameter, initial_tangent_vector, chart=None,
                 parameters=None, name=None, latex_name=None,
                 is_isomorphism=False, is_identity=False):
        r"""
        Constructs a numerical curve.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns + [x1], D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            Traceback (most recent call last):
            ...
            ValueError: Number of equations should equal codomain
             dimension.
            sage: c = M.integrated_curve(eqns, D + [x1],(t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            Traceback (most recent call last):
            ...
            ValueError: Number of velocities should equal codomain
             dimension.
            sage: c = M.integrated_curve(eqns, D, (t, -oo, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            Traceback (most recent call last):
            ...
            ValueError: Both boundaries of the interval need to be
             finite.
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), x1,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            Traceback (most recent call last):
            ...
            TypeError: x1 should be a tangent vector.
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                     name='c', parameters=[m, q, L, T])
            Traceback (most recent call last):
            ...
            TypeError: B_0 should either be a coordinate function of
             Chart (M, (x1, x2, x3)), or one of the corresponding
             velocities [Dx1, Dx2, Dx3], or the curve parameter t, or
             one of the parameters [m, q, L, T].
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:            name='c', parameters=[B_0, m, q, L, T]); c
            Integrated curve c in the 3-dimensional differentiable
             manifold M
            sage: # TestSuite(c).run() # pickling and category failed

        """

        # starting with parent class method to initialize the four last
        # arguments
        DifferentiableCurve.__init__(self, parent, name=name,
                   latex_name=latex_name, is_isomorphism=is_isomorphism,
                   is_identity=is_identity) # (coord_expression=None)

        # checking argument 'parent'
        domain = self.domain()
        if not isinstance(domain, OpenInterval):
            raise TypeError("{} is not a real interval".format(domain))
        else:
            t_min = domain.lower_bound()
            t_max = domain.upper_bound()
            if t_min == -Infinity or t_max == +Infinity:
                raise ValueError("Both boundaries of the interval " +
                                 "need to be finite.")

        # checking argument 'equations_rhs'
        codomain_dimension = self.codomain().dim()
        if len(equations_rhs) != codomain_dimension:
            raise ValueError("Number of equations should equal " +
                             "codomain dimension.")
        for eqn in equations_rhs:
            if not isinstance(eqn, Expression):
                raise TypeError("{} should be ".format(eqn) +
                                "a symbolic expression.")
        # desolve_system_rk4 called in 'solve' method is in charge of
        # raising errors about possibly remaining problems in argument
        # equations_rhs

        # checking argument 'velocities'
        if len(velocities) != codomain_dimension:
            raise ValueError("Number of velocities should equal " +
                             "codomain dimension.")

        # checking argument 'curve_parameter'
        if not isinstance(curve_parameter, Expression):
            raise TypeError("{} should be ".format(curve_parameter) +
                             "a symbolic expression.")
        # desolve_system_rk4 called in 'solve' method is in charge of
        # raising errors about possibly remaining problems in argument
        # curve_parameter

        # checking argument 'initial_tangent_vector'
        if not isinstance(initial_tangent_vector, TangentVector):
            raise TypeError("{} ".format(initial_tangent_vector) +
                            "should be a tangent vector.")

        # setting the chart
        if chart is None:
            chart = self.codomain().default_chart()

        # checking argument 'parameters'
        codomain_dimension = self.codomain().dim()
        if parameters is not None:
            for param in parameters:
                if not isinstance(param, Expression):
                    raise TypeError("{} should be ".format(param) +
                                    "a symbolic expression.")

        # checking that equations_rhs only uses the chart coordinates,
        # the velocities, the curve parameter (and the parameters if
        # there are)
        announced_variables = [coord_func for coord_func in chart[:]]
        announced_variables += [vel for vel in velocities]
        announced_variables += [curve_parameter]
        if parameters is not None:
            announced_variables += [param for param in parameters]

        equations_rhs_variables = set()
        for eqn in equations_rhs:
            equations_rhs_variables=equations_rhs_variables.union(eqn.variables())

        for rhs_variable in equations_rhs_variables:
            if rhs_variable not in announced_variables:
                str_error = "{} should either be ".format(rhs_variable)
                str_error += "a coordinate function of "
                str_error += "{}, ".format(chart)
                str_error += "or one of the corresponding velocities "
                str_error += "{}, ".format(velocities)
                str_error += "or the curve parameter "
                str_error += "{}".format(curve_parameter)
                if parameters is not None:
                    str_error += ", or one of the parameters "
                    str_error += "{}".format(parameters)
                raise TypeError(str_error + ".")

        # defining all attributes
        self._equations_rhs = list(equations_rhs) # converts to list
        # since might not already be a list (which is later required)
        self._velocities = list(velocities) # converts to list
        # since might not already be a list (which is later required)
        self._curve_parameter = curve_parameter
        self._initial_tangent_vector = initial_tangent_vector
        self._chart = chart
        self._parameters = parameters
        self._solutions = {} # dictionary containing all numerically
        # computed lists of points of the curve, the keys being chosen
        # by the user when calling method 'solve'
        self._interpolations = {} # dictionary containing lists of
        # interpolation objects, each interpolation object implementing
        # the interpolation of one of the numerical coordinate curves,
        # and the keys being chosen by the user when calling
        # method 'interpolate'

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                      parameters=[B_0, m, q, L, T]); c
            Integrated curve in the 3-dimensional differentiable
             manifold M
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:            name='c', parameters=[B_0, m, q, L, T]); c
            Integrated curve c in the 3-dimensional differentiable
             manifold M

        """

        description = "Integrated curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def system(self, verbose=True):
        r"""
        Provides a detailed description of the system defining the curve
        and returns the system defining it: chart, equations and initial
        conditions.

        INPUT:

        - ``verbose`` -- (default: ``True``) prints a detailed
          description of the curve

        OUTPUT:

        - list containing the attributes :attr:`equations_rhs`,
          :attr:`initial_tangent_vector` and :attr:`chart`

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
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
            d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
            d(x2)/dt = Dx2
            d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
            d(x3)/dt = Dx3
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
            initial_tgt_vec_comps = v0[:, initial_coord_basis] # will
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

            zip_sys = zip(chart[:],self._velocities,self._equations_rhs)
            for coord_func, velocity, eqn in zip_sys:
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)
            print(description)

        return [self._equations_rhs, v0, chart]

    def solve(self, step=0.1, method=None, solution_key=None,
              parameters_values=None, verbose=True):
        r"""
        Integrates the curve numerically over the domain of integration.

        INPUT:

        - ``step`` -- (default: ``0.1``) step of integration
        - ``method`` -- (default: ``None``) numerical scheme to use for
          the integration of the curve; algorithms available are

          * 'rk4', which uses Sage solver ``desolve_system_rk4``

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

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            sage: sol = c.solve(parameters_values={m:1, q:1, L:10, T:1})
            Traceback (most recent call last):
            ...
            ValueError: Numerical values should be provided for each of
             the parameters [B_0, m, q, L, T].
            sage: sol = c.solve(method='my method',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Traceback (most recent call last):
            ...
            ValueError: No available method of integration referred to
             as 'my method'.
            sage: sol = c.solve(
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Performing 4th order Runge-Kutta integration by default...
             Resulting list of points will be associated with the key
             'rk4' by default.
            Numerical integration completed. Resulting list of points
             was associated with the key 'rk4' (if this key already
             referred to a former numerical solution, such a solution
             was erased).
            sage: sol_mute = c.solve(verbose=False,
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            sage: sol_mute == sol
            True

        """

        if method is None:
            method = 'rk4'
            if verbose:
                print("Performing 4th order Runge-Kutta integration " +
                      "by default...")

        if solution_key is None:
            solution_key = method
            if verbose:
                print("Resulting list of points will be associated " +
                      "with the key '{}' ".format(solution_key) +
                      "by default.")

        if method == 'rk4':
            t_min = self.domain().lower_bound()
            t_max = self.domain().upper_bound()

            eqns_rhs = self._equations_rhs

            v0 = self._initial_tangent_vector
            chart = self._chart

            initial_tgt_space = v0.parent()
            initial_pt = initial_tgt_space.base_point() # retrieves
            # the initial point as the base point of the tangent space
            # to which initial tangent vector belongs
            initial_pt_coords = list(initial_pt.coordinates(chart))
            # previous line converts to list since would otherwise be a
            # tuple (yet will need to be added to [t_min] later); will
            # raise error if coordinates in chart cannot be obtained


            initial_coord_basis = chart.frame().at(initial_pt)
            initial_tgt_vec_comps = list(v0[:,initial_coord_basis])#idem

            if self._parameters is not None:
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
                    t_min=t_min.substitute(parameters_values)
                    if t_min==-Infinity or t_min==+Infinity:
                        raise ValueError("Both boundaries of the " +
                                          "interval need to be finite.")

                if isinstance(t_max, Expression):
                    t_max=t_max.substitute(parameters_values)
                    if t_max==-Infinity or t_max==+Infinity:
                        raise ValueError("Both boundaries of the " +
                                         "interval need to be finite.")

                eqns_rhs=[eqn.substitute(parameters_values) for eqn
                          in eqns_rhs]



                for i in range(len(initial_pt_coords)):
                    if isinstance(initial_pt_coords[i],Expression):
                        AUX = initial_pt_coords[i]
                        AUX = AUX.substitute(parameters_values)
                        initial_pt_coords[i] = AUX
                        # 'AUX' only used for the lines of
                        # source code to be shorter


                for i in range(len(initial_tgt_vec_comps)):
                    if isinstance(initial_tgt_vec_comps[i],Expression):
                        AUX = initial_tgt_vec_comps[i]
                        AUX = AUX.substitute(parameters_values)
                        initial_tgt_vec_comps[i] = AUX
                        # 'AUX' only used for the lines of
                        # source code to be shorter

            t_min = numerical_approx(t_min)
            t_max = numerical_approx(t_max)
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

            ics = [t_min] + initial_pt_coords + initial_tgt_vec_comps

            sol = desolve_system_rk4(self._velocities + eqns_rhs,
                                     list(chart[:]) + self._velocities,
                                     ivar = self._curve_parameter,
                                     ics = ics,
                                     end_points=[t_min, t_max],
                                     step=step)

            dim = self.codomain().dim()
            self._solutions[solution_key] = [point[0:dim+1] for point
                                             in sol]
        else:
            raise ValueError("No available method of integration " +
                             "referred to as '{}'.".format(method))

        if verbose:
            print("Numerical integration completed. " +
              "Resulting list of points was associated with the key " +
              "'{}' ".format(solution_key) +
              "(if this key already referred to a former numerical " +
              "solution, such a solution was erased).")

        return self._solutions[solution_key]

    def solution(self, solution_key=None, verbose=True):
        r"""
        Returns the solution (list of points) associated with the given
        key.

        INPUT:

        - ``solution_key`` -- (default: ``None``) key which the
          requested numerical solution is associated to ; a default
          value is chosen if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the solution returned

        OUTPUT:

        - list of the numerical points of the solution requested

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
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
            if 'rk4' in self._solutions.keys():
                solution_key = 'rk4'
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
        Interpolates the chosen numerical solution using the given
        interpolation method.

        INPUT:

        - ``solution_key`` -- (default: ``None``) key which the
          numerical solution to interpolate is associated to ; a default
          value is chosen if none is provided
        - ``method`` -- (default: ``None``) interpolation scheme to use;
          algorithms available are

          * 'cubic spline', which uses ``GSL`` via Sage class
            :class:`~sage.calculus.interpolation.Spline`

        - ``interpolation_key`` -- (default: ``None``) key which the
          resulting interpolation will be associated to ; a default
          value is given if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the interpolation in progress

        OUTPUT:

        - built interpolation object

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            sage: sol = c.solve(method='rk4', solution_key='sol_T1',
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
            if 'rk4' in self._solutions.keys():
                solution_key = 'rk4'
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
        Returns the interpolation object associated with the given key.

        INPUT:

        - ``interpolation_key`` -- (default: ``None``) key which the
          requested interpolation is associated to ; a default
          value is chosen if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the interpolation object returned

        OUTPUT:

        - requested interpolation object

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            sage: sol = c.solve(method='rk4', solution_key='sol_T1',
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
        Returns the image of the curve for the given value of the curve
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
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            sage: sol = c.solve(method='rk4', solution_key='sol_T1',
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
            [1.060743431308544, -0.2153838226258469, 1.1]
            sage: pt = c(1.1, verbose=False); pt
            [1.060743431308544, -0.2153838226258469, 1.1]

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
        Returns the vector tangent to the curve at the given curve
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

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                name='c', parameters=[B_0, m, q, L, T])
            sage: sol = c.solve(method='rk4', solution_key='sol_T1',
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
            [0.7392716344834512, -0.6734470583131389,
             0.9999999999999999]
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
    def plot(self, chart=None, ambient_coords=None, mapping=None,
             prange=None, interpolation_key=None,
             include_end_point=(True, True),
             end_point_offset=(0.001, 0.001), verbose=True, color='red',
             style='-', label_axes=True, display_tangent=False,
             color_tangent='blue', **kwds):
        r"""
        Plots the 2D or 3D projection of the curve onto the space of the
        chosen two or three ambient coordinates, based on the
        interpolation of a numerical solution previously computed.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedCurve.plot`
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
            sage: eqns = [D[1], -D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 6), v,name='c')
            sage: sol = c.solve(verbose=False)
            sage: interp = c.interpolate(verbose=False)
            sage: c_plot_2d = c.plot(ambient_coords=[x1, x2],
            ....:                 thickness=2.5,
            ....:                 display_tangent=True, plot_points=200,
            ....:                 plot_points_tangent=10, scale=0.5,
            ....:                 color='blue', color_tangent='red')
            Plotting from the interpolation associated with the key
             'cubic spline-interp-rk4' by default...
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
            eqns = [D[1], -D[0], SR(0)]
            p = M.point((0,0,0), name='p')
            Tp = M.tangent_space(p)
            v = Tp((1,0,1))
            c = M.integrated_curve(eqns, D, (t, 0, 6), v, name='c')
            c.solve(verbose=False)
            c.interpolate(verbose=False)
            c_plot_2d_1 = c.plot(ambient_coords=[x1, x2],
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
    Constructs a numerical autoparallel curve on the manifold with
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
        parameter of the curve
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in
      which the equations are given; if ``None`` the default chart
      of the manifold is assumed
    - ``parameters`` -- list of the symbolic expressions used in the
      coefficients of ``affine_connection`` other than the
      coordinates associated with the chart
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
        sage: polar.<th,ph>=S2.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: epolar = polar.frame()

    Normalizing :MATH:`e_{\phi}` provides an orthonormal basis::

        sage: ch_basis = S2.automorphism_field()
        sage: ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        sage: epolar_ON=S2.default_frame().new_frame(ch_basis,'epolar_ON')

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

    Declare the corresponding integrated autoparallel curve and display
    the differential system it satisfies::

        sage: [t, tmin, tmax] = var('t tmin tmax')
        sage: c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax),
        ....:       v, chart=polar,
        ....:       parameters=[tmin,tmax,th0,ph0,v_th0,v_ph0],name='c')
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
        d(Dth)/dt = 0
        d(ph)/dt = Dph
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
        ....:     graph2D_mercator+=c.plot(interpolation_key='interp-'+key,
        ....:                            chart=mercator, thickness=2,
        ....:                            verbose=False)

    Prepare a grid of Mercator coordinates lines, and plot the curves
    over it::

        sage: graph2D_mercator_coords=mercator.plot(chart=mercator,
        ....:                           number_values=8,color='yellow')
        sage: (graph2D_mercator + graph2D_mercator_coords).show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = S2.default_frame().new_frame(ch_basis, 'e')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                 chart=polar,parameters=[tmin,tmax,th0,ph0,v_th0,v_ph0],
                 name='c')
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
            graph2D_mercator += c.plot(interpolation_key='interp-'+key,
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
        ....:     graph3D_embedded_curves+=c.plot(interpolation_key='interp-'+key,
        ....:         mapping=euclid_embedding, thickness=5,
        ....:         display_tangent=True, scale=0.4, width_tangent=0.5,
        ....:         verbose=False)
        sage: graph3D_embedded_polar_coords = polar.plot(chart=cart,
        ....:                          mapping=euclid_embedding,
        ....:                          number_values=15, color='yellow')
        sage: graph=graph3D_embedded_curves+graph3D_embedded_polar_coords
        sage: viewer3D = 'threejs'
        sage: graph.show(viewer=viewer3D)

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = S2.default_frame().new_frame(ch_basis, 'e')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                 chart=polar,parameters=[tmin,tmax,th0,ph0,v_th0,v_ph0],
                 name='c')
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
            graph3D_embedded_curves+=c.plot(interpolation_key='interp-'+key,
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

        sage: graph2D_mercator_angle_curve=c.plot(interpolation_key='interp-angle',
        ....:         chart=mercator, thickness=1, display_tangent=True,
        ....:         scale=0.2, width_tangent=0.2, verbose=False)
        sage: graph2D_mercator_angle_curve.show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = S2.default_frame().new_frame(ch_basis, 'e')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                 chart=polar,parameters=[tmin,tmax,th0,ph0,v_th0,v_ph0],
                 name='c')
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
        graph2D_mercator_angle_curve=c.plot(interpolation_key='interp-angle',
                      chart=mercator, thickness=1, display_tangent=True,
                      scale=0.2, width_tangent=0.2, verbose=False)
        sphinx_plot(graph2D_mercator_angle_curve)

    One may eventually plot such a curve on :MATH:`\mathbb{S}^{2}`::

        sage: graph3D_embedded_angle_curve=c.plot(interpolation_key='interp-angle',
        ....:        mapping=euclid_embedding, thickness=5,
        ....:        display_tangent=True, scale=0.1, width_tangent=0.5,
        ....:        verbose=False)
        sage: graph=graph3D_embedded_angle_curve+graph3D_embedded_polar_coords
        sage: graph.show(viewer=viewer3D)

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = S2.default_frame().new_frame(ch_basis, 'e')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                 chart=polar,parameters=[tmin,tmax,th0,ph0,v_th0,v_ph0],
                 name='c')
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
        graph3D_embedded_angle_curve=c.plot(interpolation_key='interp-angle',
            mapping=euclid_embedding, thickness=5, display_tangent=True,
            scale=0.1, width_tangent=0.5, verbose=False)
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
             mapping=euclid_embedding, number_values=15, color='yellow')
        graph=graph3D_embedded_angle_curve+graph3D_embedded_polar_coords
        sphinx_plot(graph)

    """

    def __init__(self, parent, affine_connection, curve_parameter,
                 initial_tangent_vector, chart=None, parameters=None,
                 name=None, latex_name=None, is_isomorphism=False,
                 is_identity=False):
        r"""Constructs an autoparallel curve with respect to the given
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
            ....:                       name='c', parameters=[A, B]); c
            Integrated autoparallel curve c in the 3-dimensional
             differentiable manifold M

        """

        from sage.symbolic.ring import SR

        dim = parent.codomain().dim()
        i0 = parent.codomain().start_index()
        equations_rhs = []

        # setting the chart to gain access to the coordinate functions
        if chart is None:
            chart = parent.codomain().default_chart()

        coordinate_functions = chart[:]
        velocities = chart.symbolic_velocities()
        gamma = affine_connection.coef()

        for alpha in range(dim):
            rhs = SR(0)
            for mu in range(dim):
                for nu in range(dim):
                    AUX = velocities[mu] * velocities[nu]
                    rhs-= gamma[alpha+i0, mu+i0, nu+i0].expr() * AUX
                    # 'AUX' only used for the line above to be shorter
            equations_rhs += [rhs.simplify_full()]

        IntegratedCurve.__init__(self, parent, equations_rhs,
                                 velocities, curve_parameter,
                                 initial_tangent_vector, chart=chart,
                                 parameters=parameters,
                                 name=name, latex_name=latex_name,
                                 is_isomorphism=is_isomorphism,
                                 is_identity=is_identity)

        self._affine_connection = affine_connection

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

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
            ....:                                 parameters=[A, B]); c
            Integrated autoparallel curve in the 3-dimensional
             differentiable manifold M
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                       name='c', parameters=[A, B]); c
            Integrated autoparallel curve c in the 3-dimensional
             differentiable manifold M

        """

        description = "Integrated autoparallel curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def system(self, verbose=True):
        r"""
        Provides a detailed description of the system defining the
        autoparallel curve and returns the system defining it: chart,
        equations and initial conditions.

        INPUT:

        - ``verbose`` -- (default: ``True``) prints a detailed
          description of the curve

        OUTPUT:

        - list containing the attributes :attr:`equations_rhs`,
          :attr:`initial_tangent_vector` and :attr:`chart`

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
            ....:                                     parameters=[A, B])
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
            d(Dx1)/dt = -A*Dx1*Dx2*x1^2
            d(x2)/dt = Dx2
            d(Dx2)/dt = 0
            d(x3)/dt = Dx3
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
            initial_tgt_vec_comps = v0[:,initial_coord_basis] # will
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

            zip_sys = zip(chart[:],self._velocities,self._equations_rhs)
            for coord_func, velocity, eqn in zip_sys:
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)
            print(description)

        return [self._equations_rhs, v0, chart]

class IntegratedGeodesic(IntegratedAutoparallelCurve):
    r"""
    Constructs a numerical geodesic on the manifold with respect to a
    given metric.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet`
      the set of curves `\mathrm{Hom}(I, M)` to which the curve belongs
    - ``metric`` --
      :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
      metric with respect to which the curve is a geodesic
    - ``curve_parameter`` -- symbolic expression to be used as the
        parameter of the curve;
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in
      which the equations are given; if ``None`` the default chart
      of the manifold is assumed
    - ``parameters`` -- list of the symbolic expressions used in the
      coefficients of ``metric`` other than the coordinates
      associated with the chart
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
        sage: polar.<th,ph>=S2.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
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
        ....:       chart=polar,
        ....:       parameters=[tmin,tmax,th0,ph0,v_th0,v_ph0],name='c')
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
        d(Dth)/dt = Dph^2*cos(th)*sin(th)
        d(ph)/dt = Dph
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
        ....:     graph3D_embedded_geods+=c.plot(interpolation_key='interp-'+key,
        ....:                      mapping=euclid_embedding, thickness=5,
        ....:                      display_tangent=True, scale=0.3,
        ....:                      width_tangent=0.5, verbose=False)

    Plot the resulting geodesics on the grid of polar coordinates lines
    on :MATH:`\mathbb{S}^{2}` and check that these are great circles::

        sage: graph3D_embedded_polar_coords = polar.plot(chart=cart,
        ....:                          mapping=euclid_embedding,
        ....:                          number_values=15, color='yellow')
        sage: graph=graph3D_embedded_geods+graph3D_embedded_polar_coords
        sage: viewer3D = 'threejs'
        sage: graph.show(viewer=viewer3D)

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        g = S2.metric('g')
        g[1,1], g[2,2] = 1, (sin(th))**2
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar.at(p))
        c = S2.integrated_geodesic(g, (t, tmin, tmax), v, chart=polar,
                   parameters=[tmin,tmax,th0,ph0,v_th0,v_ph0], name='c')
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
            graph3D_embedded_geods+=c.plot(interpolation_key='interp-'+key,
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
                 initial_tangent_vector, chart=None, parameters=None,
                 name=None, latex_name=None, is_isomorphism=False,
                 is_identity=False):

        r"""Constructs a geodesic curve with respect to the given metric
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
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v, name='c',
            ....:                            parameters=[A]); c
            Integrated geodesic c in the 2-dimensional differentiable
             manifold S^2

        """

        IntegratedAutoparallelCurve.__init__(self, parent,
                                   metric.connection(), curve_parameter,
                                   initial_tangent_vector, chart=chart,
                                   parameters=parameters, name=name,
                                   latex_name=latex_name,
                                   is_isomorphism=is_isomorphism,
                                   is_identity=is_identity)

        self._metric = metric

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

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
            ....:                            parameters=[A]); c
            Integrated geodesic in the 2-dimensional differentiable
             manifold S^2
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v, name='c',
            ....:                            parameters=[A]); c
            Integrated geodesic c in the 2-dimensional differentiable
             manifold S^2

        """

        description = "Integrated geodesic "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def system(self, verbose=True):
        r"""
        Returns the system defining the geodesic : chart, equations and
        initial conditions

        INPUT:

        - ``verbose`` -- (default: ``True``) prints a detailed
          description of the curve

        OUTPUT:

        - list containing the attributes :attr:`equations_rhs`,
          :attr:`initial_tangent_vector` and :attr:`chart`

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
            ....:                            parameters=[A])
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
            d(Dtheta)/dt = Dphi^2*cos(theta)*sin(theta)
            d(phi)/dt = Dphi
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
            initial_tgt_vec_comps=v0[:,initial_coord_basis]#will
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

            zip_sys = zip(chart[:],self._velocities,self._equations_rhs)
            for coord_func, velocity, eqn in zip_sys:
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)
            print(description)

        return [self._equations_rhs, v0, chart]
