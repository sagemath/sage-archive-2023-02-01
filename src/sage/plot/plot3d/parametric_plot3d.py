"""
3d Parametric Plots
"""

from parametric_surface import ParametricSurface
from shapes2 import line3d
from texture import Texture

def parametric_plot3d(f, urange, vrange=None, plot_points="automatic", **kwds):
    """
    Return a parametric three-dimensional space curve or surface.

    INPUTS:
    There are four ways to call this function:
        * parametric_plot3d([f_x, f_y, f_z], (u_min, u_max)):
          f_x, f_y, f_z are three functions and u_min and u_max
          are real numbers

        * parametric_plot3d([f_x, f_y, f_z], (u, u_min, u_max))
          f_x, f_y, f_z can be viewed as functions of u

        * parametric_plot3d([f_x, f_y, f_z], (u_min, u_max), (v_min, v_max))
          f_x, f_y, f_z are each functions of two variables

        * parametric_plot3d([f_x, f_y, f_z], (u, u_min, u_max), (v, v_min, v_max))
          f_x, f_y, f_z can be viewed as functions of u and v

    The INPUTS are as follows:
        f -- a 3-tuple of functions or expressions
        urange -- a 2-tuple (u_min, u_max) or a 3-tuple (u, u_min, u_max)
        vrange -- (optional -- only used for surfaces) a 2-tuple (u_min, u_max)
                  or a 3-tuple (u, u_min, u_max)
        plot_points -- (default: "automatic", which is 75 for curves and
                       [15,15] for surfaces) initial number of sample
                       points in each parameter; an integer for a curve,
                       and a pair of integers for a surface.

    NOTES:
      * By default for a curve any points where f_x, f_y, or f_z do
        not evaluate to a real number are skipped.
      * Currently for a surface f_x, f_y, and f_z have to be defined
        everywhere. This will change.

    EXAMPLES:
    We demonstrate each of the four ways to call this function.

        1. A space curve defined by three functions of 1 variable:
            sage: parametric_plot3d( (sin, cos, lambda u: u/10), (0, 20))

        Note above the lambda function, which creates a callable Python function
        that sends u to u/10.

        2. Next we draw the same plot as above, but using symbolic functions:
            sage: u = var('u')
            sage: parametric_plot3d( (sin(u), cos(u), u/10), (u, 0, 20))

        3. We draw a parametric surface using 3 Python functions (defined using
           lambda):
            sage: f = (lambda u,v: cos(u), lambda u,v: sin(u)+cos(v), lambda u,v: sin(v))
            sage: parametric_plot3d(f, (0, 2*pi), (-pi, pi))

        4. The same surface, but where the defining functions are symbolic:
            sage: u, v = var('u,v')
            sage: parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)), (u, 0, 2*pi), (v, -pi, pi))

        We increase the number of plot points, and make the surface green and transparent:
            sage: parametric_plot3d((cos(u), sin(u) + cos(v), sin(v)), (u, 0, 2*pi), (v, -pi, pi), color='green', opacity=0.1, plot_points=[30,30])

    We call the space curve function but with polynomials instead of
    symbolic variables.

        sage: R.<t> = RDF[]
        sage: parametric_plot3d( (t, t^2, t^3), (t, 0, 3) )

    Next we plot the same curve, but because we use (0, 3) instead of (t, 0, 3),
    each polynomial is viewed as a callable function of one variable:

        sage: parametric_plot3d( (t, t^2, t^3), (0, 3) )

    We do a plot but mix a symbolic input, and an integer:
        sage: t = var('t')
        sage: parametric_plot3d( (1, sin(t), cos(t)), (t, 0, 3) )

    We plot two interlinked tori:
        sage: u, v = var('u,v')
        sage: f1 = (4+(3+cos(v))*sin(u), 4+(3+cos(v))*cos(u), 4+sin(v))
        sage: f2 = (8+(3+cos(v))*cos(u), 3+sin(v), 4+(3+cos(v))*sin(u))
        sage: p1 = parametric_plot3d(f1, (u,0,2*pi), (v,0,2*pi), texture="red")
        sage: p2 = parametric_plot3d(f2, (u,0,2*pi), (v,0,2*pi), texture="blue")
        sage: p1 + p2

    A Conchoid:
        sage: u, v = var('u,v')
        sage: k = 1.2; k_2 = 1.2; a = 1.5
        sage: f = (k^u*(1+cos(v))*cos(u), k^u*(1+cos(v))*sin(u), k^u*sin(v)-a*k_2^u)
        sage: parametric_plot3d(f, (u,0,6*pi), (v,0,2*pi), plot_points=[40,40], texture=(0,0.5,0))
    """
    # TODO:
    #   * Surface -- behavior of functions not defined everywhere -- see note above
    #   * Iterative refinement


    # boundary_style -- (default: None) how boundary lines are drawn for a surface
    # color_function -- (default: "automatic") how to determine the color of curves and surfaces
    # color_function_scaling -- (default: True) whether to scale the input to color_function
    # exclusions -- (default: "automatic") u points or (u,v) conditions to exclude.
    #         (E.g., exclusions could be a function e = lambda u, v: False if u < v else True
    # exclusions_style -- (default: None) what to draw at excluded points
    # max_recursion -- (default: "automatic") maximum number of recursive subdivisions,
    #                   when ...
    # mesh -- (default: "automatic") how many mesh divisions in each direction to draw
    # mesh_functions -- (default: "automatic") how to determine the placement of mesh divisions
    # mesh_shading -- (default: None) how to shade regions between mesh divisions
    # plot_range -- (default: "automatic") range of values to include

    if isinstance(f, (list,tuple)) and len(f) > 0 and isinstance(f[0], (list,tuple)):
        return sum([parametric_plot3d(v, urange, vrange, plot_points, **kwds) for v in f])

    if not isinstance(f, (tuple, list)) or len(f) != 3:
        raise ValueError, "f must be a list or tuple of length 3"

    if vrange is None:
        if plot_points == "automatic":
            plot_points = 75
        G = parametric_plot3d_curve(f, urange, plot_points, **kwds)
    else:
        if plot_points == "automatic":
            plot_points = [15,15]
        G = parametric_plot3d_surface(f, urange, vrange, plot_points, **kwds)
    G._set_extra_kwds(kwds)
    return G

def parametric_plot3d_curve(f, urange, plot_points, **kwds):
    from sage.plot.plot import var_and_list_of_values
    plot_points = int(plot_points)
    u, vals = var_and_list_of_values(urange, plot_points)
    w = []
    fail = 0

    if u is None:
        f_x, f_y, f_z = f
        for t in vals:
            try:
                w.append((float(f_x(t)), float(f_y(t)), float(f_z(t))))
            except TypeError:
                fail += 1
    else:
        f_x, f_y, f_z = [ensure_subs(m) for m in f]
        for t in vals:
            try:
                w.append((float(f_x.subs({u:t})), float(f_y.subs({u:t})),
                          float(f_z.subs({u:t}))))
            except TypeError:
                fail += 1

    if fail > 0:
        print "WARNING: Failed to evaluate parametric plot at %s points"%fail
    return line3d(w, **kwds)

def parametric_plot3d_surface(f, urange, vrange, plot_points, **kwds):
    if not isinstance(plot_points, (list, tuple)) or len(plot_points) != 2:
        raise ValueError, "plot_points must be a tuple of length 2"
    points0, points1 = plot_points

    from sage.plot.plot import var_and_list_of_values
    u, u_vals = var_and_list_of_values(urange, int(points0))
    v, v_vals = var_and_list_of_values(vrange, int(points1))

    if u is None:
        if not v is None:
            raise ValueError, "both ranges must specify a variable or neither must"
        # nothing to do
        f_x, f_y, f_z = f
    else:
        if v is None:
            raise ValueError, "both ranges must specify a variable or neither must"
        f0, f1, f2 = [ensure_subs(w) for w in f]
        def f_x(uu,vv):
            return float(f0.subs({u:uu, v:vv}))
        def f_y(uu,vv):
            return float(f1.subs({u:uu, v:vv}))
        def f_z(uu,vv):
            return float(f2.subs({u:uu, v:vv}))

    def g(x,y):
        # Change to use fast callable float symbolic expressions later
        return (float(f_x(x,y)), float(f_y(x,y)), float(f_z(x,y)))

    return ParametricSurface(g, (u_vals, v_vals), **kwds)

def ensure_subs(f):
    if not hasattr(f, 'subs'):
        from sage.calculus.all import SR
        return SR(f)
    return f
