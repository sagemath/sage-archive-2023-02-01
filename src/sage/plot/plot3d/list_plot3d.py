"""
List Plots
"""

from sage.matrix.all import is_Matrix, matrix
from sage.rings.all  import RDF

def list_plot3d(v, interpolation_type='default', texture="automatic", point_list=None,**kwds):
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

    OPTIONAL KEYWORDS:

        interpolation_type - 'linear', 'nn' (nearest neighbor), 'spline'

             'linear'  will perform linear interpolation

             The option 'nn' will interpolate by averaging the value
             of the nearest neighbors, this produces an interpolating function that is smoother than a linear interpolation, it
             has one derivative everywhere except at the sample points.

             The option 'spline' interpolates using a bivariate B-spline.

             When v is a matrix the default is to use linear interpolation,
             when v is a list of points the default is nearest neighbor.

        degree - an integer between 1 and 5, controls the degree of spline
             used for spline interpolation. For data that is highly oscillatory
             use higher values

        point_list - If point_list=True is passed, then if the array is a
             list of lists of length three, it will be treated  as an array of
             points rather than a 3xn array.

        num_points - Number of points to sample interpolating function in each direction. By default
                     for an nxn array this is n.


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

   Though it doesn't change the shap of the graph, increasing num_points can increase the clarity of the
   graph
        sage: list_plot3d(m, texture='yellow', frame_aspect_ratio=[1,1,1/3],num_points=40)

    We can change the interpolation type
        sage: list_plot3d(m, texture='yellow', interpolation_type='nn',frame_aspect_ratio=[1,1,1/3])

    We can make this look better by increasing the number of samples
        sage: list_plot3d(m, texture='yellow', interpolation_type='nn',frame_aspect_ratio=[1,1,1/3],num_points=40)

    Lets try a spline
        sage: list_plot3d(m, texture='yellow', interpolation_type='spline',frame_aspect_ratio=[1,1,1/3])

    That spline doesn't capture the oscillation very well, lets try a higher degree spline
        sage: list_plot3d(m, texture='yellow', interpolation_type='spline', degree=5, frame_aspect_ratio=[1,1,1/3])


    We plot a list of lists:
        sage: show(list_plot3d([[1, 1, 1, 1], [1, 2, 1, 2], [1, 1, 3, 1], [1, 2, 1, 4]]))

    We plot a list of points:
        As a first example  we can extract the (x,y,z) coordinates from the above example and
        make a list plot out of it. By default we do linear interpolation.

        sage: l=[]
        sage: for i in range(6):
        ...      for j in range(6):
        ...         l.append((float(i*pi/5),float(j*pi/5),m[i,j]))
        sage: list_plot3d(l,texture='yellow')


        Note that the points do not have to be regularly sampled. For example

        sage: l=[]
        sage: for i in range(-5,5):
        ...    for j in range(-5,5):
        ...      l.append((normalvariate(0,1),normalvariate(0,1),normalvariate(0,1)))
        sage: list_plot3d(l,interpolation_type='nn',texture='yellow',num_points=100)


    """
    import numpy
    if texture == "automatic":
        texture = "lightblue"
    if is_Matrix(v):
        if interpolation_type=='default' or interpolation_type=='linear' and not kwds.has_key('num_points'):
            return list_plot3d_matrix(v, texture=texture,  **kwds)
        else:
            l=[]
            for i in xrange(v.nrows()):
                for j in xrange(v.ncols()):
                    l.append((i,j,v[i,j]))
            return list_plot3d_tuples(l,interpolation_type,texture,**kwds)

    if type(v)==numpy.ndarray:
        return list_plot3d(matrix(v),interpolation_type,texture,**kwds)

    if isinstance(v, list):
        if len(v) == 0:
            # return empty 3d graphic
            from base import Graphics3d
            return Graphics3d()
        elif isinstance(v[0],tuple) or point_list==True and len(v[0]) == 3:
            return list_plot3d_tuples(v,interpolation_type,texture=texture, **kwds)
        else:
            return list_plot3d_array_of_arrays(v, interpolation_type,texture, **kwds)
    raise TypeError, "v must be a matrix or list"

def list_plot3d_matrix(m, texture, **kwds):
    from parametric_surface import ParametricSurface
    f = lambda i,j: (i,j,float(m[int(i),int(j)]))
    G = ParametricSurface(f, (range(m.nrows()), range(m.ncols())), texture=texture, **kwds)
    G._set_extra_kwds(kwds)
    return G

def list_plot3d_array_of_arrays(v, interpolation_type,texture, **kwds):
    m = matrix(RDF, len(v), len(v[0]), v)
    G = list_plot3d(m,interpolation_type,texture, **kwds)
    G._set_extra_kwds(kwds)
    return G

def list_plot3d_tuples(v,interpolation_type, texture, **kwds):
    import delaunay
    import numpy
    import scipy
    from random import random
    from scipy import interpolate
    from scipy import stats
    from plot3d import plot3d

    x=[float(p[0]) for p in v]
    y=[float(p[1]) for p in v]
    z=[float(p[2]) for p in v]

    corr_matrix=stats.corrcoef(x,y)

    if corr_matrix[0,1] > .9 or corr_matrix[0,1] <  -.9:
        # If the x,y coordinates lie in a one-dimensional subspace
        # The scipy delauney code segfaults
        # We compute the correlation of the x and y coordinates
        # and add small random noise to avoid the problem
        # if it looks like there is an issue

        ep=float(.000001)
        x=[float(p[0])+random()*ep for p in v]
        y=[float(p[1])+random()*ep for p in v]


    xmin=float(min(x))
    xmax=float(max(x))
    ymin=float(min(y))
    ymax=float(max(y))

    num_points= kwds['num_points'] if kwds.has_key('num_points') else int(4*numpy.sqrt(len(x)))
                                          #arbitrary choice - assuming more or less a nxn grid of points
                                          # x should have n^2 entries. We sample 4 times that many points.



    if interpolation_type == 'linear':

        T= delaunay.Triangulation(x,y)
        f=T.linear_interpolator(z)
        f.default_value=0.0
        j=numpy.complex(0,1)
        vals=f[ymin:ymax:j*num_points,xmin:xmax:j*num_points]
        from parametric_surface import ParametricSurface

        def g(x,y):
            i=round( (x-xmin)/(xmax-xmin)*(num_points-1) )
            j=round( (y-ymin)/(ymax-ymin)*(num_points-1) )
            z=vals[int(j),int(i)]
            return (x,y,z)


        G = ParametricSurface(g, (list(numpy.r_[xmin:xmax:num_points*j]), list(numpy.r_[ymin:ymax:num_points*j])), texture=texture, **kwds)
        G._set_extra_kwds(kwds)
        return G


    if interpolation_type == 'nn'  or interpolation_type =='default':

        T=delaunay.Triangulation(x,y)
        f=T.nn_interpolator(z)
        f.default_value=0.0
        j=numpy.complex(0,1)
        vals=f[ymin:ymax:j*num_points,xmin:xmax:j*num_points]
        from parametric_surface import ParametricSurface
        def g(x,y):
            i=round( (x-xmin)/(xmax-xmin)*(num_points-1) )
            j=round( (y-ymin)/(ymax-ymin)*(num_points-1) )
            z=vals[int(j),int(i)]
            return (x,y,z)

        G = ParametricSurface(g, (list(numpy.r_[xmin:xmax:num_points*j]), list(numpy.r_[ymin:ymax:num_points*j])), texture=texture, **kwds)
        G._set_extra_kwds(kwds)
        return G

    if interpolation_type =='spline':
        from plot3d import plot3d
        kx=kwds['kx'] if kwds.has_key('kx') else 3
        ky=kwds['ky'] if kwds.has_key('ky') else 3
        if kwds.has_key('degree'):
            kx=kwds['degree']
            ky=kwds['degree']

        s=kwds['smoothing'] if kwds.has_key('smoothing') else len(x)-numpy.sqrt(2*len(x))
        s=interpolate.bisplrep(x,y,z,[int(1)]*len(x),xmin,xmax,ymin,ymax,kx=kx,ky=ky,s=s)
        f=lambda x,y: interpolate.bisplev(x,y,s)
        return plot3d(f,(xmin,xmax),(ymin,ymax),texture=texture,plot_points=[num_points,num_points],**kwds)
