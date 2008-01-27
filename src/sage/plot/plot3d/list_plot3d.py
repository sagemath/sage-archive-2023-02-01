"""
List Plots
"""

from sage.matrix.all import is_Matrix, matrix
from sage.rings.all  import RDF

def list_plot3d(v, texture="automatic", **kwds):
    """
    A 3-dimensional plot of a surface defined by the list $v$ of
    points in 3-dimensional space.

    INPUT:
        v -- something that defines a set of points in 3 space,
             for example:
                 * a matrix
                 * a list of 3-tuples
                 * a list of lists (all of the same length) -- this
                   is treated the same as a matrix.
        texture -- (default: "automatic"), solid light blue
        **kwds -- all other arguments are passed to the surface function

    OUTPUT:
        a 3d plot

    EXAMPLES:
    We plot a matrix that illustrates summation modulo $n$.
        sage: n = 5; list_plot3d(matrix(RDF,n,[(i+j)%n for i in [1..n] for j in [1..n]]))

    We plot a matrix of values of sin.
        sage: pi = float(pi)
        sage: m = matrix(RDF, 6, [sin(i^2 + j^2) for i in [0,pi/5,..,pi] for j in [0,pi/5,..,pi]])
        sage: list_plot3d(m, texture='yellow', frame_aspect_ratio=[1,1,1/3])

    We plot a list of lists:
        sage: show(list_plot3d([[1, 1, 1, 1], [1, 2, 1, 2], [1, 1, 3, 1], [1, 2, 1, 4]]))

    We plot a list of points:
        NOT IMPLEMENTED YET

    """
    if texture == "automatic":
        texture = "lightblue"
    if is_Matrix(v):
        return list_plot3d_matrix(v, texture=texture,  **kwds)
    if isinstance(v, list):
        if len(v) == 0:
            # return empty 3d graphic
            from base import Graphics3d
            return Graphics3d()
        elif isinstance(v[0],tuple) and len(v[0]) == 3:
            return list_plot3d_tuples(v, texture=texture, **kwds)
        else:
            return list_plot3d_array_of_arrays(v, texture=texture, **kwds)
    raise TypeError, "v must be a matrix or list"

def list_plot3d_matrix(m, texture, **kwds):
    from parametric_surface import ParametricSurface
    f = lambda i,j: (i,j,float(m[int(i),int(j)]))
    G = ParametricSurface(f, (range(m.nrows()), range(m.ncols())), texture=texture, **kwds)
    G._set_extra_kwds(kwds)
    return G

def list_plot3d_array_of_arrays(v, texture, **kwds):
    m = matrix(RDF, len(v), len(v[0]), v)
    G = list_plot3d_matrix(m, texture=texture, **kwds)
    G._set_extra_kwds(kwds)
    return G

def list_plot3d_tuples(v, texture, **kwds):
    raise NotImplementedError
