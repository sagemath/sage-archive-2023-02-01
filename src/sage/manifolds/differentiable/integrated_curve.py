r"""
Integrated Curves in Manifolds

Given a differentiable manifold `M`, an *integrated curve* curve in `M` is a
differentiable curve constructed as a numerical solution to a system of
second order differential equations.

Integrated curves are implemented by :class:`IntegratedCurve`, which the
classes :class:`IntegratedAutoparallelCurve` and :class:`IntegratedGeodesic`
inherit.

"""

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
    :MATH:`t \mapsto (x_{1}(t), ..., x_{n}(t))` constructed as a numerical
    solution to a system of second order differential equations satisfied by
    the coordinate curves :MATH:`t \mapsto x_{i}(t)`.
    
    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet`
      the set of curves `\mathrm{Hom}(I, M)` to which the curve belongs
    - ``equations_rhs`` -- list of the right-hand sides of the equations on
      the velocities only (the term *velocity* referring to the derivatives
      :MATH:`d x_{i} / dt` of the coordinate curves)
    - ``velocities`` -- list of the symbolic expressions used in ``equations_rhs``
      to denote the velocities
    - ``curve_parameter`` -- symbolic expression used in ``equations_rhs``
      to denote the parameter of the curve (denoted :MATH:`t` in the
      descriptions above)
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      the curve; if none is provided, ``name`` will be used
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

        \mathbf{B}(t,\mathbf{x}) = \frac{B_{0}t}{T} \exp \left( - \frac{ x_{1}^{2} + x_{2}^{2} }{ L^{2} }  \right) \mathbf{e_{3}}.
    
    Equations of motion are:
    
    .. MATH::
    
        \ddot{x}_{1}(t) &= \frac{qB(t, \mathbf{x}(t))}{m} \dot{x}_{2}(t)    \\        
        \ddot{x}_{2}(t) &= - \frac{qB(t, \mathbf{x}(t))}{m} \dot{x}_{1}(t)   \\        
        \ddot{x}_{3}(t) &= 0
    
    Start with declaring a chart on a 3-dimensional manifold and the various
    symbolic expressions denoting the equations, velocities and parameters::    

        sage: M = Manifold(3, 'M')
        sage: X.<x1,x2,x3> = M.chart()
        sage: var('t B_0 m q L T')
        (t, B_0, m, q, L, T)
        sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
        sage: D = X.symbolic_velocities() ; D
        [Dx1, Dx2, Dx3]
        sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
        
    Set the initial conditions::
            
        sage: p = M.point((0,0,0), name='p')
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((1,0,1))
        
    Declare an integrated curve and display information relative to it::
    
        sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
        ....:                        parameters=[B_0, m, q, L, T]) ; c
        Integrated curve c in the 3-dimensional differentiable manifold M
        sage: c.system()
        Curve c in the 3-dimensional differentiable manifold M integrated over
         the Real interval (0, 5) as a solution to the following system, written
         with respect to the Chart (M, (x1, x2, x3)):

         Initial point: Point p on the 3-dimensional differentiable manifold M
         with coordinates [0, 0, 0] in Chart (M, (x1, x2, x3))
         Initial tangent vector: Tangent vector at Point p on the 3-dimensional
         differentiable manifold M with components [1, 0, 1] in Chart (M, (x1, x2, x3))

         d(x1)/dt = Dx1
         d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
         d(x2)/dt = Dx2
         d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
         d(x3)/dt = Dx3
         d(Dx3)/dt = 0

    Generate a solution of the system and an interpolation of this solution::
    
        sage: c.solve(step=0.2, parameters_values={B_0:1, m:1, q:1, L:10, T:1},
        ....:         solution_key='carac time 1')
        Performing 4th order Runge-Kutta integration by default...
         Numerical integration completed. Resulting list of points was associated
         with the key 'carac time 1' (if this key already referred to a former
         numerical solution, such a solution was erased).
        sage: c.interpolate(solution_key='carac time 1', interpolation_key='interp 1')
        Performing cubic spline interpolation by default...
         Interpolation completed and associated with the key 'interp 1'
         (if this key already referred to a former interpolation,
         such an interpolation was erased).
         
    Such an interpolation is required to evaluate the curve and the vector
    tangent to the curve for any value of the curve parameter::

        sage: c(1.9)
        Evaluating point coordinates from the interpolation associated with
         the key 'interp 1' by default...
         [1.3776707219621374, -0.9000776970132945, 1.9]
        sage: v = c.tangent_vector_eval_at(4.3)
        Evaluating tangent vector components from the interpolation associated
         with the key 'interp 1' by default...
        sage: v
        Tangent vector at Point on the 3-dimensional differentiable manifold M
        sage: v.display()
        -0.9303968397216424 d/dx1 - 0.3408080563014475 d/dx2 + 1.0000000000000004 d/dx3

    Plotting a numerical solution (with or without its tangent vector field)
    also requires the solution to be interpolated at least once::

        sage: c_plot_2d_1 = c.plot(ambient_coords=[x1, x2], interpolation_key='interp 1',
        ....:                      thickness=2.5, display_tangent=True, plot_points=200, 
        ....:                      plot_points_tangent=10, scale=0.5, color='blue',
        ....:                      color_tangent='red')
        sage: c_plot_2d_1.show()

    .. PLOT::
    
        M = Manifold(3, 'M')
        X = M.chart('x1 x2 x3')
        var('x1 x2 x3 t B_0 m q L T')
        B = B_0*t/T*exp(-(x1**2 + x2**2)/L**2)
        D = X.symbolic_velocities() ; D
        eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
                               parameters=[B_0, m, q, L, T])
        c.solve(step=0.2, parameters_values={B_0:1, m:1, q:1, L:10, T:1},
                solution_key='carac time 1')
        c.interpolate(solution_key='carac time 1', interpolation_key='interp 1')
        c_plot_2d_1 = c.plot(ambient_coords=[x1, x2], interpolation_key='interp 1',
                             thickness=2.5, display_tangent=True, plot_points=200,
                             plot_points_tangent=10, scale=0.5, color='blue',
                             color_tangent='red')
        sphinx_plot(c_plot_2d_1)
        
    An instance of :class:`IntegratedCurve` may store several numerical solutions
    and interpolations::

        sage: c.solve(step=0.2, parameters_values={B_0:1, m:1, q:1, L:10, T:100},
        ....:         solution_key='carac time 100')
        Performing 4th order Runge-Kutta integration by default...
         Numerical integration completed. Resulting list of points was associated
         with the key 'carac time 100' (if this key already referred to a former
         numerical solution, such a solution was erased).
        sage: c.interpolate(solution_key='carac time 100', interpolation_key='interp 100')
        Performing cubic spline interpolation by default...
         Interpolation completed and associated with the key 'interp 100'
         (if this key already referred to a former interpolation, such an
         interpolation was erased).
        sage: c_plot_3d_100 = c.plot(interpolation_key='interp 100', thickness=2.5,
        ....:                        display_tangent=True, plot_points=200,
        ....:                        plot_points_tangent=10, scale=0.5,
        ....:                        color='green', color_tangent='orange')
        sage: c_plot_3d_1 = c.plot(interpolation_key='interp 1', thickness=2.5,
        ....:                      display_tangent=True, plot_points=200,
        ....:                      plot_points_tangent=10, scale=0.5,
        ....:                      color='blue', color_tangent='red')
        sage: (c_plot_3d_1 + c_plot_3d_100).show()

    .. PLOT::
    
        M = Manifold(3, 'M')
        X = M.chart('x1 x2 x3')
        var('x1 x2 x3 t B_0 m q L T')
        B = B_0*t/T*exp(-(x1**2 + x2**2)/L**2)
        D = X.symbolic_velocities() ; D
        eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
                               parameters=[B_0, m, q, L, T])
        c.solve(step=0.2, parameters_values={B_0:1, m:1, q:1, L:10, T:1},
                solution_key='carac time 1')
        c.interpolate(solution_key='carac time 1', interpolation_key='interp 1')
        c.solve(step=0.2, parameters_values={B_0:1, m:1, q:1, L:10, T:100},
                solution_key='carac time 100')
        c.interpolate(solution_key='carac time 100', interpolation_key='interp 100') 
        c_plot_3d_1 = c.plot(interpolation_key='interp 1', thickness=2.5,
                             display_tangent=True, plot_points=200,
                             plot_points_tangent=10, scale=0.5, color='blue',
                             color_tangent='red')
        c_plot_3d_100 = c.plot(interpolation_key='interp 100', thickness=2.5,
                               display_tangent=True, plot_points=200,
                               plot_points_tangent=10, scale=0.5, color='green',
                               color_tangent='orange')
        sphinx_plot(c_plot_3d_1 + c_plot_3d_100)        
    
    """
    
    def __init__(self, parent, equations_rhs, velocities, curve_parameter,
                 initial_tangent_vector, chart=None, parameters=None, name=None,
                 latex_name=None, is_isomorphism=False, is_identity=False):
                 
        r"""
        Constructs a numerical curve.
        
        TESTS::
        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns + [x1], D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            ValueError: Number of equations should equal codomain dimension.
            sage: c = M.integrated_curve(eqns, D + [x1], (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            ValueError: Number of velocities should equal codomain dimension.
            sage: c = M.integrated_curve(eqns, D, (t, -oo, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            ValueError: Both boundaries of the interval need to be finite.
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), x1, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            TypeError: x1 should be a tangent vector.
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
                                         parameters=[m, q, L, T])
            TypeError: B_0 should either be a coordinate function of B_0, or
             one the corresponding velocities [Dx1, Dx2, Dx3], or the curve
             parameter t, or one of the parameters [m, q, L, T].            
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T]) ; c
            Integrated curve c in the 3-dimensional differentiable manifold M
            sage: # TestSuite(c).run() # pickling and category failed

        """        
        
        # starting with parent class method to initialize the four last arguments
        DifferentiableCurve.__init__(self, parent, name=name, latex_name=latex_name,
                                     is_isomorphism=is_isomorphism,
                                     is_identity=is_identity) # (coord_expression=None)
        
        # checking argument 'parent'
        domain = self.domain()
        if not isinstance(domain, OpenInterval):
            raise TypeError("{} is not a real interval".format(domain))
        else:
            t_min = domain.lower_bound()
            t_max = domain.upper_bound()
            if t_min == -Infinity or t_max == +Infinity:
                raise ValueError("Both boundaries of the interval need to be finite.")

        # checking argument 'equations_rhs'       
        codomain_dimension = self.codomain().dim()
        if len(equations_rhs) != codomain_dimension:
            raise ValueError("Number of equations should equal codomain dimension.")
        for eqn in equations_rhs:
            if not isinstance(eqn, Expression):
                raise TypeError("{} should be a symbolic expression.".format(eqn))
        # desolve_system_rk4 called in 'solve' method is in charge of raising
        # errors about possibly remaining problems in argument equations_rhs

        # checking argument 'velocities'        
        if len(velocities) != codomain_dimension:
            raise ValueError("Number of velocities should equal codomain dimension.")

        # checking argument 'curve_parameter'
        if not isinstance(curve_parameter, Expression):
            raise TypeError("{} should be a symbolic expression.".format(curve_parameter))
        # desolve_system_rk4 called in 'solve' method is in charge of raising
        # errors about possibly remaining problems in argument curve_parameter
        
        # checking argument 'initial_tangent_vector'
        if not isinstance(initial_tangent_vector, TangentVector):
            raise TypeError("{} should be a tangent vector.".format(initial_tangent_vector))

        # setting the chart
        if chart is None:
            chart = self.codomain().default_chart()
            
        # checking argument 'parameters'       
        codomain_dimension = self.codomain().dim()
        if parameters is not None:
            for param in parameters:
                if not isinstance(param, Expression):
                    raise TypeError("{} should be a symbolic expression.".format(param))
            
        # checking that equations_rhs only uses the chart coordinates, the velocities,
        # the curve parameter and the parameters if there are
        announced_variables = [coord_func for coord_func in chart[:]] 
        announced_variables += [vel for vel in velocities] + [curve_parameter]
        if parameters is not None:
            announced_variables += [param for param in parameters]

        equations_rhs_variables = set()
        for eqn in equations_rhs:
            equations_rhs_variables = equations_rhs_variables.union(eqn.variables())

        for rhs_variable in equations_rhs_variables:
            if rhs_variable not in announced_variables:
                str_error = "{} should either be ".format(rhs_variable)
                str_error += "a coordinate function of {}, ".format(rhs_variable, chart)
                str_error += "or one the corresponding velocities {}, ".format(velocities)
                str_error += "or the curve parameter {}".format(curve_parameter)
                if parameters is not None:
                    str_error += ", or one of the parameters {}".format(parameters)                    
                raise TypeError(str_error + ".")
                                
        # defining all attributes
        self._equations_rhs = [eqn for eqn in equations_rhs] # converts to list
        # since might not already be a list (which is later required)
        self._velocities = [vel for vel in velocities] # converts to list since
        # might not already be a list (which is later required)
        self._curve_parameter = curve_parameter
        self._initial_tangent_vector = initial_tangent_vector
        self._chart = chart
        self._parameters = parameters                      
        self._solutions = {} # dictionary containing all numerically computed
        # lists of points of the curve, the keys being chosen by the user when
        # calling method 'solve'
        self._interpolations = {} # dictionary containing lists of interpolation
        # objects, each interpolation object implementing the interpolation of
        # one of the numerical coordinate curves, and the keys being chosen by
        # the user when calling method 'interpolate'



    def _repr_(self):
        r"""        
        Returns a string representation of ``self``.

        TESTS::
        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                        parameters=[B_0, m, q, L, T]) ; c
            Integrated curve in the 3-dimensional differentiable manifold M
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c'
            ....:                        parameters=[B_0, m, q, L, T]) ; c        
            Integrated curve c in the 3-dimensional differentiable manifold M
        
        """

        description = "Integrated curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description



    def system(self):
        r"""
        Provides a detailed description of the system defining the curve.
        
        TESTS::
        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            sage: c.system()
            Curve c in the 3-dimensional differentiable manifold M integrated
             over the Real interval (0, 5) as a solution to the following system,
             written with respect to the Chart (M, (x1, x2, x3)):

             Initial point: Point p on the 3-dimensional differentiable manifold M
             with coordinates [0, 0, 0] in Chart (M, (x1, x2, x3))
             Initial tangent vector: Tangent vector at Point p on the 3-dimensional
             differentiable manifold M with components [1, 0, 1] in Chart (M, (x1, x2, x3))

             d(x1)/dt = Dx1
             d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
             d(x2)/dt = Dx2
             d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
             d(x3)/dt = Dx3
             d(Dx3)/dt = 0
            
        """
        
        initial_tangent_space = self._initial_tangent_vector.parent()
        initial_point = initial_tangent_space.base_point() # gets the initial point
        # as the base point of the tangent space to which initial tangent vector belongs
        initial_point_coordinates = initial_point.coordinates(self._chart) # will
        # raise error if coordinates in chart are not known        
        initial_point_coordinates = [coord for coord in initial_point_coordinates] # converts
        # to list since was previously a tuple
        
        initial_coordinate_basis = self._chart.frame().at(initial_point)
        initial_tangent_vector_components = self._initial_tangent_vector[:, initial_coordinate_basis] # will
        # raise error if components in coordinate basis are not known        

        description = "Curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {} integrated over the {} ".format(self.codomain(),
                                                                  self.domain())
                                                                  
        description += "as a solution to the following system, written with "
        description += "respect to the {}:\n\n".format(self._chart)        
        description += "Initial point: {} with coordinates {} in {}\n".format(initial_point,
                                                                    initial_point_coordinates,
                                                                    self._chart)
        description += "Initial tangent vector: {} with components {} in {}\n\n".format(self._initial_tangent_vector,
                                                                     initial_tangent_vector_components,
                                                                     self._chart)
        
        for coord_func, velocity, eqn in zip(self._chart[:], self._velocities, self._equations_rhs):
            description += "d({})/d{} = {}\n".format(coord_func, self._curve_parameter, velocity)
            description += "d({})/d{} = {}\n".format(velocity, self._curve_parameter, eqn)        
        print(description)
        
        return [self._equations_rhs, self._initial_tangent_vector, self._chart]
        
        

    def solve(self, step=0.1, method=None, solution_key=None, parameters_values=None):
        r"""
        Integrates the curve over the domain of integration.
        
        Integration scheme 'rk4' uses Sage solver *desolve_system_rk4*.

        TESTS::
        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])       
            sage: c.solve(parameters_values={m:1, q:1, L:10, T:1})
            ValueError: Numerical values should be provided for each of the
             parameters [B_0, m, q, L, T].
            sage: c.solve(method='my method', parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            ValueError: No available method of integration referred to as 'my method'.
            sage: c.solve(parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Performing 4th order Runge-Kutta integration by default...
             Resulting list of points will be associated with the key 'rk4' by
             default.
             Numerical integration completed. Resulting list of points was
             associated with the key 'rk4' (if this key already referred to a
             former numerical solution, such a solution was erased).

        """

        if method is None:
            method = 'rk4'
            print('Performing 4th order Runge-Kutta integration by default...')
            
        if solution_key is None:
            solution_key = method
            print("Resulting list of points will be associated with the key '{}' by default.".format(solution_key))
        
        if method == 'rk4':
            if self._parameters is not None:
                if parameters_values is None or len(parameters_values) != len(self._parameters):
                    raise ValueError("Numerical values should be provided for each of the parameters {}.".format(self._parameters))
                eqns_rhs = [eqn.substitute(parameters_values) for eqn in self._equations_rhs] # will
                # raise error if any element of parameters_values is not numerical
            else:
                eqns_rhs = self._equations_rhs
           
            coordinate_functions_list = [coord_func for coord_func in self._chart[:]]

            t_min = numerical_approx(self.domain().lower_bound())
            t_max = numerical_approx(self.domain().upper_bound())

            initial_tangent_space = self._initial_tangent_vector.parent()
            initial_point = initial_tangent_space.base_point() # gets the initial
            # point as the base point of the tangent space to which the initial tangent vector belongs
            initial_point_coordinates = initial_point.coordinates(self._chart) # will
            # raise error if coordinates in chart cannot get known        
            initial_point_coordinates = [numerical_approx(coord) for coord
                                         in initial_point_coordinates] # converts
            # to list since was previously a tuple (so it can be added to [t_min] below),
            # and gets numerical value in case some coordinates contained expressions such as pi
            
            initial_coordinate_basis = self._chart.frame().at(initial_point)
            initial_tangent_vector_components = self._initial_tangent_vector[:, initial_coordinate_basis] # will
            # raise error if components in coordinate basis cannot get known
            initial_tangent_vector_components = [numerical_approx(comp) for comp
                                                 in initial_tangent_vector_components] # should
            # already be a list, but the components need to be completely numerical as well
            
            sol = desolve_system_rk4(self._velocities + eqns_rhs,
                                     coordinate_functions_list + self._velocities,
                                     ivar = self._curve_parameter,
                                     ics=[t_min] + initial_point_coordinates + initial_tangent_vector_components,
                                     end_points=[t_min, t_max], step=step)
            
            dim = self.codomain().dim()
            self._solutions[solution_key] = [point[0:dim+1] for point in sol]
        else:
            raise ValueError("No available method of integration referred to as '{}'.".format(method))

        print("Numerical integration completed. " +
              "Resulting list of points was associated with the key '{}' ".format(solution_key) +
              "(if this key already referred to a former numerical solution, " +
              "such a solution was erased).")
        
        return self._solutions[solution_key]



    def solution(self, solution_key=None):
        r"""
        Returns the list of points associated with the given key.
        
        TESTS::        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            sage: c.solve(parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Performing 4th order Runge-Kutta integration by default...
             Resulting list of points will be associated with the key 'rk4' by
             default.
             Numerical integration completed. Resulting list of points was
             associated with the key 'rk4' (if this key already referred to a
             former numerical solution, such a solution was erased).
            sage: c.solution(solution_key='my solution')
            ValueError: No existing key 'my solution' referring to any numerical
             solution.
            sage: c.solution()
            Returning the numerical solution associated with the key 'rk4' by
             default...

        """
        
        if solution_key is None:
            if 'rk4' in self._solutions.keys():
                solution_key = 'rk4'
            else:
                solution_key = self._solutions.keys()[0] # will raise error if
                # self._solutions empty
            print("Returning the numerical solution associated with the key '{}' by default...".format(solution_key))
        elif solution_key not in self._solutions.keys():
            raise ValueError("No existing key '{}' referring to any numerical solution.".format(solution_key))
            
        return self._solutions[solution_key]



    def interpolate(self, solution_key=None, method=None, interpolation_key=None):
        r"""        
        Interpolates the chosen numerical solution using the given interpolation
        method. 
               
        Interpolation scheme 'cubic spline' uses *GSL* via Sage class 
        :class:`~sage.calculus.interpolation.Spline`.
        
        TESTS::        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            sage: c.solve(method='rk4', solution_key='sol1',
            ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Numerical integration completed. Resulting list of points was associated
             with the key 'sol1' (if this key already referred to a former
             numerical solution, such a solution was erased).
            sage: c.interpolate(solution_key='my solution')
            ValueError: No existing key 'my solution' referring to any numerical
             solution.
            sage: c.interpolate(solution_key='sol1', method='my method')
            ValueError: No available method of interpolation referred to as 'my method'.
            sage: c.interpolate(solution_key='sol1')
            Performing cubic spline interpolation by default...
             Resulting interpolation will be associated with the key
             'cubic spline-interp-sol1' by default.
             Interpolation completed and associated with the key
             'cubic spline-interp-sol1' (if this key already referred to a
             former interpolation, such an interpolation was erased).
            sage: c.interpolate(method='cubic spline', solution_key='sol1',
                                interpolation_key='interp1')
            Interpolation completed and associated with the key 'interp1' (if
             this key already referred to a former interpolation, such an
             interpolation was erased).
 
        """
        
        if solution_key is None:
            if 'rk4' in self._solutions.keys():
                solution_key = 'rk4'
            else:
                solution_key = self._solutions.keys()[0] # will raise error if
                # self._solutions empty
            print("Interpolating the numerical solution associated with the key '{}' by default...".format(solution_key))
        elif solution_key not in self._solutions.keys():
            raise ValueError("No existing key '{}' referring to any numerical solution.".format(solution_key))
        
        if method is None:
            method = 'cubic spline'
            print("Performing cubic spline interpolation by default...") 

        if interpolation_key is None:
            interpolation_key = "{}-interp-{}".format(method, solution_key)
            print("Resulting interpolation will be associated with the key '{}' by default.".format(interpolation_key))
        
        if method=='cubic spline':
            self._interpolations[interpolation_key] = []
            dim = self.codomain().dim()
            for i in range(dim):
                coordinate_curve = []
                for point in self._solutions[solution_key]:
                    coordinate_curve += [[point[0], point[i+1]]]
                self._interpolations[interpolation_key] += [Spline(coordinate_curve)]
        else:
            raise ValueError("No available method of interpolation referred to as '" +
                             method + "'.")

        print("Interpolation completed and associated with the key '{}' ".format(interpolation_key) +
              "(if this key already referred to a former interpolation, " +
              "such an interpolation was erased).")
        
        return self._interpolations[interpolation_key]


       
    def interpolation(self, interpolation_key=None):
        r"""
        Returns the interpolation object associated with the given key.
        
        TESTS::        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            sage: c.solve(method='rk4', solution_key='sol1',
            ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Numerical integration completed. Resulting list of points was associated
             with the key 'sol1' (if this key already referred to a former
             numerical solution, such a solution was erased).
            sage: c.interpolate(method='cubic spline', solution_key='sol1',
                                interpolation_key='interp1')
            Interpolation completed and associated with the key 'interp1' (if
             this key already referred to a former interpolation, such an
             interpolation was erased).
            sage: c.interpolation(interpolation_key='my interp')
            ValueError: No existing key 'my interp' referring to any interpolation.
            sage: c.interpolation()
            Returning the interpolation associated with the key 'interp1' by
             default...

        """

        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0] # will
                # raise error if self._interpolations empty
            print("Returning the interpolation associated with the key '{}' ".format(interpolation_key)
                  + "by default...")
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key '{}' referring to any interpolation.".format(interpolation_key)) 
        
        return self._interpolations[interpolation_key]



    def __call__(self, curve_parameter_value, interpolation_key=None):
        r"""
        Returns the image of the curve for the given value of the curve
        parameter, using the chosen interpolation.        
        
        TESTS::        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            sage: c.solve(method='rk4', solution_key='sol1',
            ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Numerical integration completed. Resulting list of points was associated
             with the key 'sol1' (if this key already referred to a former
             numerical solution, such a solution was erased).
            sage: c.interpolate(method='cubic spline', solution_key='sol1',
                                interpolation_key='interp1')
            Interpolation completed and associated with the key 'interp1' (if
             this key already referred to a former interpolation, such an
             interpolation was erased).
            sage: c(1.1, interpolation_key='my interp')
            ValueError: No existing key 'my interp' referring to any
             interpolation.
            sage: c(1.1)
            Evaluating point coordinates from the interpolation associated
             with the key 'interp1' by default...
             [1.060743431308544, -0.2153838226258469, 1.1]
        
        """
        
        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0] # will
                # raise error if self._interpolations empty
            print("Evaluating point coordinates from the interpolation " +
                  "associated with the key '{}' ".format(interpolation_key) +
                  "by default...")
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key '{}' referring to any interpolation.".format(interpolation_key))            
        
        interpolation = self._interpolations[interpolation_key]
        
        if isinstance(interpolation[0], Spline): #partial test, in case
        # future interpolation objects do not contain lists of instances of
        # the Spline class
            interpolated_coordinates = [coordinate_curve_spline(curve_parameter_value)
                                        for coordinate_curve_spline in interpolation]
            return interpolated_coordinates

        raise TypeError("Unexpected type of interpolation object.")
    
    
    
    def tangent_vector_eval_at(self, curve_parameter_value, interpolation_key=None):
        r"""        
        Returns the vector tangent to the curve at the given curve parameter
        with components evaluated from the given interpolation.
        
        TESTS::        
            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t B_0 m q L T')
            (t, B_0, m, q, L, T)
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
            ....:                        parameters=[B_0, m, q, L, T])
            sage: c.solve(method='rk4', solution_key='sol1',
            ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Numerical integration completed. Resulting list of points was associated
             with the key 'sol1' (if this key already referred to a former
             numerical solution, such a solution was erased).
            sage: c.interpolate(method='cubic spline', solution_key='sol1',
                                interpolation_key='interp1')
            Interpolation completed and associated with the key 'interp1' (if
             this key already referred to a former interpolation, such an
             interpolation was erased).
        
        """

        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0] # will
                # raise error if self._interpolations empty
            print("Evaluating tangent vector components from the interpolation " +
                  "associated with the key '{}' ".format(interpolation_key) +
                  "by default...")
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key '{}' referring to any interpolation.".format(interpolation_key))
        
        interpolation = self._interpolations[interpolation_key]
        
        if isinstance(interpolation[0], Spline): #partial test, in case
        # future interpolation objects do not contain lists of instances of
        # the Spline class
            interpolated_coordinates = [coordinate_curve_spline(curve_parameter_value)
                                              for coordinate_curve_spline in interpolation]
            M = self.codomain()
            p = M.point(interpolated_coordinates, chart=self._chart)
            Tp = M.tangent_space(p)
                                                          
            evaluated_tangent_vector_components = [coordinate_curve_spline.derivative(curve_parameter_value, order=1)
                                                   for coordinate_curve_spline in interpolation]
            basis = self._chart.frame().at(p)
            v = Tp(evaluated_tangent_vector_components, basis=basis)
            return v

        raise TypeError("Unexpected type of interpolation object.")



    @options(thickness=1, plot_points=75, plot_points_tangent=10, aspect_ratio='automatic', scale=1)
    def plot(self, display_tangent=False, color_tangent='blue', 
             interpolation_key=None, ambient_coords=None, prange=None,
             include_end_point=(True, True), end_point_offset=(0.001, 0.001),
             color='red', style='-', label_axes=True, **kwds):
        r"""
        Plots the 2D or 3D projection of the curve onto the space of the chosen
        two or three ambient coordinates, based on the interpolation of a
        numerical solution previously computed."
        
        TESTS::        
            sage: #TO DO
        
        """

        #
        # Get the @options from kwds
        #
        thickness = kwds.pop('thickness')
        plot_points = kwds.pop('plot_points')
        aspect_ratio = kwds.pop('aspect_ratio')
        
        #
        # Coordinates of the chart w.r.t. which the curve is plotted
        #
        if ambient_coords is None:
            ambient_coords = self._chart[:]  # all chart coordinates are used
        n_pc = len(ambient_coords)
        if n_pc != 2 and n_pc !=3:
            raise ValueError("the number of coordinates involved in the " +
                             "plot must be either 2 or 3, not {}".format(n_pc))
        # indices of plot coordinates
        ind_pc = [self._chart[:].index(pc) for pc in ambient_coords]
        
        #
        # Parameter range for the plot
        #        
        param_min = numerical_approx(self.domain().lower_bound())
        param_max = numerical_approx(self.domain().upper_bound())
                
        if prange is None:
            prange = (param_min, param_max)
        elif not isinstance(prange, (tuple, list)):
            raise TypeError("{} is neither a tuple nor a list".format(prange))
        elif len(prange) != 2:
            raise ValueError("the argument prange must be a tuple/list " +
                             "of 2 elements")

        tmin = numerical_approx(prange[0])
        tmax = numerical_approx(prange[1])

        if tmin < param_min or tmin > param_max or tmax < param_min or tmax > param_max:
            raise ValueError("Parameter range should be a subinterval of the curve domain ({}).".format(self.domain()))
                    
        if not include_end_point[0]:
            tmin = tmin + end_point_offset[0]

        if not include_end_point[1]:
            tmax = tmax - end_point_offset[1]
        
        tmin = numerical_approx(tmin)
        tmax = numerical_approx(tmax)        

        #
        # List of points for the plot curve
        #
        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0] # will raise error if
                # self._interpolations empty
            print("Plotting from the interpolation associated " +
                  "with the key '{}' ".format(interpolation_key) +
                  "by default...")
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key '{}' referring to any interpolation.".format(interpolation_key))            
        
        interpolation = self._interpolations[interpolation_key]
        
        if isinstance(interpolation[0], Spline): #partial test, in case
        # future interpolation objects do not contain lists of instances of
        # the Spline class
            plot_curve = []
            dt = (tmax - tmin) / (plot_points - 1)
            t = tmin
            

            for i in range(plot_points):
                plot_curve.append( [interpolation[j](t) for j in ind_pc] )
                t += dt                   
                                         
            if display_tangent:
                from sage.plot.graphics import Graphics
                from sage.plot.arrow import arrow2d
                from sage.plot.plot3d.shapes import arrow3d
                
                scale = kwds.pop('scale')
                plot_points_tangent = kwds.pop('plot_points_tangent', None)
                if plot_points_tangent is None:
                    plot_points_tangent = plot_points
                    print("Plotting as many tangent vectors as points.")                
                
                plot_vectors = Graphics()
                dt = (tmax - tmin) / (plot_points_tangent - 1)
                t = tmin
                                
                for i in range(plot_points_tangent):
                    # interpolated ambient coordinates:
                    xp = [interpolation[j](t) for j in ind_pc]
                                        
                    # tangent vector ambiant components evaluated from the interpolation:
                    vcomp = [coordinate_curve_spline.derivative(t, order=1)
                             for coordinate_curve_spline in interpolation]                    
                    
                    coord_tail = xp
                    coord_head = [xp[i] + scale*vcomp[i] for i in ind_pc]
                    
                    if coord_head != coord_tail:
                        if n_pc == 2:
                            plot_vectors += arrow2d(tailpoint=coord_tail,
                                                    headpoint=coord_head,
                                                    color=color_tangent)
                        else:
                            plot_vectors += arrow3d(coord_tail, coord_head,
                                                    color=color_tangent)
                    t += dt
                return plot_vectors + DifferentiableCurve._graphics(self, plot_curve, 
                                             ambient_coords, thickness=thickness,
                                             aspect_ratio=aspect_ratio, color=color,
                                             style=style, label_axes=label_axes)
                                             
                                             
            return DifferentiableCurve._graphics(self, plot_curve, ambient_coords, thickness=thickness,
                                         aspect_ratio=aspect_ratio, color=color,
                                         style=style, label_axes=label_axes) 
        
        raise TypeError("Unexpected type of interpolation object.")        



class IntegratedAutoparallelCurve(IntegratedCurve):
    r"""    
    Constructs a numerical autoparallel curve.
        
    INPUT:
                TO DO
            
    """

    def __init__(self, parent, affine_connection, curve_parameter, initial_tangent_vector,
                 chart=None, parameters=None, name=None, latex_name=None,
                 is_isomorphism=False, is_identity=False):
                 
        r"""Construct the autoparallel curve with respect to the given affine connection with the given initial tangent vector."""
        
        from sage.symbolic.ring import SR

        dim = parent.codomain().dim()
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
                    rhs -= gamma[alpha, mu, nu].expr()*velocities[mu]*velocities[nu]
            equations_rhs += [rhs.simplify_full()]
        
        IntegratedCurve.__init__(self, parent, equations_rhs, velocities, curve_parameter,
                                 initial_tangent_vector, chart=chart, parameters=parameters,
                                 name=name, latex_name=latex_name, is_isomorphism=is_isomorphism,
                                  is_identity=is_identity)
        
        self._affine_connection = affine_connection


    def _repr_(self):
        r"""        
        Returns a string representation of ``self``.

        TESTS::        
            sage: #TO DO

        """

        description = "Integrated autoparallel curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description   


        
    def system(self):
        r"""
        Returns the system defining the autoparallel curve : chart, equations and initial conditions
        
        TESTS::        
            sage: #TO DO        
        
        """
        
        initial_tangent_space = self._initial_tangent_vector.parent()
        initial_point = initial_tangent_space.base_point() # gets the initial
        # point as the base point of the tangent space which initial tangent vector belongs to
        initial_point_coordinates = initial_point.coordinates(self._chart) # will
        # raise error if coordinates in chart are not known        
        initial_point_coordinates = [coord for coord in initial_point_coordinates] # converts
        # to list since was previously a tuple
        
        initial_coordinate_basis = self._chart.frame().at(initial_point)
        initial_tangent_vector_components = self._initial_tangent_vector[:, initial_coordinate_basis] # will
        # raise error if components in coordinate basis are not known        

        description = "Autoparallel curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {} equipped with the {}, and integrated over the {} ".format(self.codomain(), self._affine_connection, self.domain())        
        description += "as a solution to the following equations, written with respect to the {}:\n\n".format(self._chart)
        
        description += "Initial point: {} with coordinates {} in {}\n".format(initial_point,
                                                                             initial_point_coordinates, self._chart)
        description += "Initial tangent vector: {} with components {} in {}\n\n".format(self._initial_tangent_vector,
                                                                                              initial_tangent_vector_components, self._chart)
        
        for coord_func, velocity, eqn in zip(self._chart[:], self._velocities, self._equations_rhs):
            description += "d({})/d{} = {}\n".format(coord_func, self._curve_parameter, velocity)
            description += "d({})/d{} = {}\n".format(velocity, self._curve_parameter, eqn)        
        print(description)
        
        return [self._equations_rhs, self._initial_tangent_vector, self._chart]



class IntegratedGeodesic(IntegratedAutoparallelCurve):
    r"""    
    Constructs a numerical geodesic on the manifold.
        
    INPUT:
                TO DO
            
    """

    def __init__(self, parent, metric, curve_parameter, initial_tangent_vector,
                 chart=None, parameters=None, name=None, latex_name=None,
                 is_isomorphism=False, is_identity=False):
                 
        r"""Construct the geodesic curve with respect to the given metric with the given initial tangent vector."""
        
        
        IntegratedAutoparallelCurve.__init__(self, parent, metric.connection(), curve_parameter,
                                 initial_tangent_vector, chart=chart, parameters=parameters,
                                 name=name, latex_name=latex_name, is_isomorphism=is_isomorphism,
                                  is_identity=is_identity)                                 
        
        self._metric = metric



    def _repr_(self):
        r"""        
        Returns a string representation of ``self``.

        TESTS::        
            sage: #TO DO

        """

        description = "Integrated geodesic "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description   


        
    def system(self):
        r"""
        Returns the system defining the geodesic : chart, equations and initial conditions
        
        TESTS::        
            sage: #TO DO
        
        """
        
        initial_tangent_space = self._initial_tangent_vector.parent()
        initial_point = initial_tangent_space.base_point() # gets the initial
        # point as the base point of the tangent space which initial tangent vector belongs to
        initial_point_coordinates = initial_point.coordinates(self._chart) # will
        # raise error if coordinates in chart are not known        
        initial_point_coordinates = [coord for coord in initial_point_coordinates] # converts
        # to list since was previously a tuple
        
        initial_coordinate_basis = self._chart.frame().at(initial_point)
        initial_tangent_vector_components = self._initial_tangent_vector[:, initial_coordinate_basis] # will
        # raise error if components in coordinate basis are not known        

        description = "Geodesic "
        if self._name is not None:
            description += self._name + " "
        description += "in the {} equipped with the {}, and integrated over the {} ".format(self.codomain(), self._metric, self.domain())        
        description += "as a solution to the following geodesic equations, written with respect to the {}:\n\n".format(self._chart)
        
        description += "Initial point: {} with coordinates {} in {}\n".format(initial_point,
                                                                             initial_point_coordinates, self._chart)
        description += "Initial tangent vector: {} with components {} in {}\n\n".format(self._initial_tangent_vector,
                                                                                              initial_tangent_vector_components, self._chart)
        
        for coord_func, velocity, eqn in zip(self._chart[:], self._velocities, self._equations_rhs):
            description += "d({})/d{} = {}\n".format(coord_func, self._curve_parameter, velocity)
            description += "d({})/d{} = {}\n".format(velocity, self._curve_parameter, eqn)        
        print(description)
        
        return [self._equations_rhs, self._initial_tangent_vector, self._chart]
