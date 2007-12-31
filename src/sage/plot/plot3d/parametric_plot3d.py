"""
3d Parametric Plots
"""

from parametric_surface import ParametricSurface
from shapes2 import line3d

def parametric_plot3d(f, range0, range1=None, plot_points=50, **kwds):
    if range1 is None:
        return parametric_plot3d_curve(f, range0, plot_points, **kwds)
    else:
        return parametric_plot3d_curve(f, range0, range1, plot_points, **kwds)

def parametric_plot3d_curve(f, range0, plot_points, **kwds):
    v = list_of_values(range0, plot_points)
    w = [f(t) for t in v]
    return line3d(w, **kwds)

def parametric_plot3d_surface(f, range0, range1, plot_points, **kwds):
    if isinstance(plot_points, (list, tuple)):
        if len(plot_points) != 2:
            raise ValueError, "plot_points must be an integer or list of 2 integers"
        points0, points1 = plot_points
    else:
        points0 = plot_points; point1 = plot_points

    range0 = list_of_values(range0, points=int(points0))
    range1 = list_of_values(range1, points=int(points1))

    if isinstance(f, (tuple, list)):
        if len(f) != 3:
            raise ValueError, "f must be either a function of two variables or a 3-tuple or list of length 3."
        def g(x,y):
            return (f[0](x,y), f[1](x,y), f[2](x,y))
    else:
        # deal with callable symbolic expressions, etc.
        g = f

    return ParametricSurface(g, (range0, range1), **kwds)



def list_of_values(v, points):
    if points < 1: points = 1
    if isinstance(v, list):
        return [float(t) for t in v]
    elif isinstance(v, tuple):
        if len(v) == 2:
            return float_range(v[0], v[1], float(v[1]-v[0])/points)
        elif len(v) == 3:
            return float_range(v[0], v[1], float(v[1]-v[0])/v[2])
    raise TypeError, "parametric value range must be a list of 2 or 3-tuple."

def float_range(a, b, step):
    a = float(a); b = float(b); step = float(step)
    v = [a]
    w = a + step
    while w <= b:
        v.append(w)
        w += step
    if w < b:
        v.append(b)
    return v
