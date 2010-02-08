r"""
Plotting Functions.

EXAMPLES::

    sage: def f(x,y):
    ...       return math.sin(y*y+x*x)/math.sqrt(x*x+y*y+.0001)
    ...
    sage: P = plot3d(f,(-3,3),(-3,3), adaptive=True, color=rainbow(60, 'rgbtuple'), max_bend=.1, max_depth=15)
    sage: P.show()

::

    sage: def f(x,y):
    ...       return math.exp(x/5)*math.sin(y)
    ...
    sage: P = plot3d(f,(-5,5),(-5,5), adaptive=True, color=['red','yellow'])
    sage: from sage.plot.plot3d.plot3d import axes
    sage: S = P + axes(6, color='black')
    sage: S.show()

We plot "cape man"::

    sage: S = sphere(size=.5, color='yellow')

::

    sage: from sage.plot.plot3d.shapes import Cone
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)

::

    sage: S += sphere((.45,-.1,.15), size=.1, color='white') + sphere((.51,-.1,.17), size=.05, color='black')
    sage: S += sphere((.45, .1,.15),size=.1, color='white') + sphere((.51, .1,.17), size=.05, color='black')
    sage: S += sphere((.5,0,-.2),size=.1, color='yellow')
    sage: def f(x,y): return math.exp(x/5)*math.cos(y)
    sage: P = plot3d(f,(-5,5),(-5,5), adaptive=True, color=['red','yellow'], max_depth=10)
    sage: cape_man = P.scale(.2) + S.translate(1,0,0)
    sage: cape_man.show(aspect_ratio=[1,1,1])

AUTHORS:

- Tom Boothby: adaptive refinement triangles

- Josh Kantor: adaptive refinement triangles

- Robert Bradshaw (2007-08): initial version of this file

- William Stein (2007-12, 2008-01): improving 3d plotting

- Oscar Lazo, William Cauchois (2009-2010): Adding coordinate transformations
"""


#TODO:
#    -- smooth triangles

#*****************************************************************************
#      Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from tri_plot import TrianglePlot
from index_face_set import IndexFaceSet
from shapes import arrow3d
from base import Graphics3dGroup
from sage.plot.colors import rainbow
from texture import Texture, is_Texture

from sage.ext.fast_eval import fast_float_arg, fast_float

from sage.functions.trig import cos, sin

class _CoordTrans(object):
    """
    This abstract class encapsulates a new coordinate system for plotting.
    Sub-classes must implement the ``gen_transform`` method which, given
    symbolic variables to use, generates a 3-tuple of functions in terms of
    those variables that can be used to find the cartesian (X, Y, and Z)
    coordinates for any point in this space.
    """

    def __init__(self, indep_var, dep_vars):
        """
        INPUT:

         - ``indep_var`` - The independent variable (the function value will be
           substituted for this).

         - ``dep_vars`` - A list of dependent variables (the parameters will be
           substituted for these).

        TESTS::

            sage: from sage.plot.plot3d.plot3d import _CoordTrans
            sage: _CoordTrans('y', ['x'])
            Unknown coordinate system (y in terms of x)
        """
        if hasattr(self, 'all_vars'):
            if set(self.all_vars) != set(dep_vars + [indep_var]):
                raise ValueError, 'not all variables were specified for ' + \
                                  'this coordinate system'
        self.indep_var = indep_var
        self.dep_vars = dep_vars

    def gen_transform(self, **kwds):
        """
        Generate the transformation for this coordinate system in terms of the
        specified variables (which should be keywords).
        """
        raise NotImplementedError

    def to_cartesian(self, func, params):
        """
        Returns a 3-tuple of functions, parameterized over ``params``, that
        represents the cartesian coordinates of the value of ``func``.

        INPUT:

         - ``func`` - A function in this coordinate space. Corresponds to the
           independent variable.

         - ``params`` - The parameters of func. Correspond to the dependent
           variables.

        EXAMPLE::

            sage: from sage.plot.plot3d.plot3d import _ArbCoordTrans
            sage: x, y, z = var('x y z')
            sage: T = _ArbCoordTrans((x + y, x - y, z), z)
            sage: f(x, y) = x * y
            sage: T.to_cartesian(f, [x, y])
            [x + y, x - y, x*y]
        """

        from sage.symbolic.expression import is_Expression
        from sage.calculus.calculus import is_RealNumber
        from sage.rings.integer import is_Integer
        if any([is_Expression(func), is_RealNumber(func), is_Integer(func)]):
            return self.gen_transform(**{
                self.indep_var: func,
                self.dep_vars[0]: params[0],
                self.dep_vars[1]: params[1]
            })
        else:
            # func might be a lambda or a Python callable; this makes it slightly
            # more complex.
            import sage.symbolic.ring
            indep_var_dummy = sage.symbolic.ring.var(self.indep_var)
            dep_var_dummies = sage.symbolic.ring.var(','.join(self.dep_vars))
            transformation = self.gen_transform(**{
                self.indep_var: indep_var_dummy,
                self.dep_vars[0]: dep_var_dummies[0],
                self.dep_vars[1]: dep_var_dummies[1]
            })
            def subs_func(t):
                return lambda x,y: t.subs({
                    indep_var_dummy: func(x, y),
                    dep_var_dummies[0]: x,
                    dep_var_dummies[1]: y
                })
            return map(subs_func, transformation)

    _name = 'Unknown coordinate system'

    def __repr__(self):
        return '%s (%s in terms of %s)' % \
          (self._name, self.indep_var, ', '.join(self.dep_vars))

class _ArbCoordTrans(_CoordTrans):
    """
    An arbitrary coordinate system transformation.
    """

    def __init__(self, custom_trans, fvar):
        """
        INPUT:

         - ``custom_trans`` - A 3-tuple of transformations. This will be returned
           almost unchanged by ``gen_transform``, except the function variable
           will be substituted.

         - ``fvar`` - The function variable.
        """
        super(_ArbCoordTrans, self).__init__('f', ['u', 'v'])
        self.custom_trans = custom_trans
        self.fvar = fvar

    def gen_transform(self, f=None, u=None, v=None):
        """
        EXAMPLE::

            sage: from sage.plot.plot3d.plot3d import _ArbCoordTrans
            sage: x, y, z = var('x y z')
            sage: T = _ArbCoordTrans((x + y, x - y, z), x)

        The independent and dependent variables don't really matter in the case
        of an arbitrary transformation (since it is already in terms of its own
        variables), so default values are provided::

            sage: T.indep_var
            'f'
            sage: T.dep_vars
            ['u', 'v']

        Finally, an example of gen_transform()::

            sage: T.gen_transform(f=z)
            [y + z, -y + z, z]
        """
        # We don't need to do anything with u and v here because self.custom_trans
        # should already be in terms of those variables.
        return [t.subs({self.fvar: f}) for t in self.custom_trans]

class Spherical(_CoordTrans):
    """
    A spherical coordinate system for use with ``plot3d(transformation=...)``
    where the position of a point is specified by three numbers:

     - the *radial distance* (``r``),

     - the *elevation angle* (``theta``),

     - and the *azimuth angle* (``phi``).

    These three variables must be specified in the constructor.

    EXAMPLES:

    Construct a spherical transformation for a function ``r`` in terms of
    ``theta`` and ``phi``::

        sage: T = Spherical('r', ['phi', 'theta'])

    If we construct some concrete variables, we can get a transformation::

        sage: r, phi, theta = var('r phi theta')
        sage: T.gen_transform(r=r, theta=theta, phi=phi)
        (r*sin(theta)*cos(phi), r*sin(phi)*sin(theta), r*cos(theta))

    Use with plot3d on a made-up function::

        sage: plot3d(phi * theta, (phi, 0, 1), (theta, 0, 1), transformation=T)

    To graph a function ``theta`` in terms of ``r`` and ``phi``, you would use::

        sage: Spherical('theta', ['r', 'phi'])
        Spherical coordinate system (theta in terms of r, phi)

    See also ``spherical_plot3d`` for more examples of plotting in spherical
    coordinates.
    """

    all_vars = ['r', 'theta', 'phi']

    _name = 'Spherical coordinate system'

    def gen_transform(self, r=None, theta=None, phi=None):
        """
        EXAMPLE::

            sage: T = Spherical('r', ['theta', 'phi'])
            sage: T.gen_transform(r=var('r'), theta=var('theta'), phi=var('phi'))
            (r*sin(theta)*cos(phi), r*sin(phi)*sin(theta), r*cos(theta))
        """
        return (r * sin(theta) * cos(phi),
                r * sin(theta) * sin(phi),
                r * cos(theta))

class Cylindrical(_CoordTrans):
    """
    A cylindrical coordinate system for use with ``plot3d(transformation=...)``
    where the position of a point is specified by three numbers:

     - the *radial distance* (``rho``),

     - the *angular position* or *azimuth* (``phi``),

     - and the *height* or *altitude* (``z``).

    These three variables must be specified in the constructor.

    EXAMPLES:

    Construct a cylindrical transformation for a function ``rho`` in terms of
    ``phi`` and ``z``::

        sage: T = Cylindrical('rho', ['phi', 'z'])

    If we construct some concrete variables, we can get a transformation::

        sage: rho, phi, z = var('rho phi z')
        sage: T.gen_transform(rho=rho, phi=phi, z=z)
        (rho*cos(phi), rho*sin(phi), z)

    Use with plot3d on a made-up function::

        sage: plot3d(phi * z, (phi, 0, 1), (z, 0, 1), transformation=T)

    To graph a function ``z`` in terms of ``phi`` and ``rho`` you would use::

        sage: Cylindrical('z', ['phi', 'rho'])
        Cylindrical coordinate system (z in terms of phi, rho)

    See also ``cylindrical_plot3d`` for more examples of plotting in cylindrical
    coordinates.
    """

    _name = 'Cylindrical coordinate system'

    all_vars = ['rho', 'phi', 'z']

    def gen_transform(self, rho=None, phi=None, z=None):
        """
        EXAMPLE::

            sage: T = Cylindrical('z', ['phi', 'rho'])
            sage: T.gen_transform(rho=var('rho'), phi=var('phi'), z=var('z'))
            (rho*cos(phi), rho*sin(phi), z)
        """
        return (rho * cos(phi),
                rho * sin(phi),
                z)

class TrivialTriangleFactory:
    def triangle(self, a, b, c, color = None):
        return [a,b,c]
    def smooth_triangle(self, a, b, c, da, db, dc, color = None):
        return [a,b,c]

import parametric_plot3d
def plot3d(f, urange, vrange, adaptive=False, transformation=None, **kwds):
    """
    INPUT:


    -  ``f`` - a symbolic expression or function of 2
       variables

    -  ``urange`` - a 2-tuple (u_min, u_max) or a 3-tuple
       (u, u_min, u_max)

    -  ``vrange`` - a 2-tuple (v_min, v_max) or a 3-tuple
       (v, v_min, v_max)

    -  ``adaptive`` - (default: False) whether to use
       adaptive refinement to draw the plot (slower, but may look better).
       This option does NOT work in conjuction with a transformation
       (see below).

    -  ``mesh`` - bool (default: False) whether to display
       mesh grid lines

    -  ``dots`` - bool (default: False) whether to display
       dots at mesh grid points

    -  ``plot_points`` - (default: "automatic") initial number of sample
       points in each direction; an integer or a pair of integers


    - ``transformation`` - (default: None) a transformation to apply. May be a
      4-tuple (x_func, y_func, z_func, fvar) where the first 3 items indicate a
      transformation to cartesian coordinates (from your coordinate system) in
      terms of u, v, and the function variable fvar (for which the value of f
      will be substituted). May also be a predefined coordinate system
      transformation like Spherical or Cylindrical.

    .. note::

       ``mesh`` and ``dots`` are not supported when using the Tachyon
       raytracer renderer.

    EXAMPLES: We plot a 3d function defined as a Python function::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2))

    We plot the same 3d function but using adaptive refinement::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2), adaptive=True)

    Adaptive refinement but with more points::

        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2), adaptive=True, initial_depth=5)

    We plot some 3d symbolic functions::

        sage: var('x,y')
        (x, y)
        sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
        sage: plot3d(sin(x*y), (x, -pi, pi), (y, -pi, pi))

    We give a plot with extra sample points::

        sage: var('x,y')
        (x, y)
        sage: plot3d(sin(x^2+y^2),(x,-5,5),(y,-5,5), plot_points=200)
        sage: plot3d(sin(x^2+y^2),(x,-5,5),(y,-5,5), plot_points=[10,100])

    A 3d plot with a mesh::

        sage: var('x,y')
        (x, y)
        sage: plot3d(sin(x-y)*y*cos(x),(x,-3,3),(y,-3,3), mesh=True)

    Two wobby translucent planes::

        sage: x,y = var('x,y')
        sage: P = plot3d(x+y+sin(x*y), (x,-10,10),(y,-10,10), opacity=0.87, color='blue')
        sage: Q = plot3d(x-2*y-cos(x*y),(x,-10,10),(y,-10,10),opacity=0.3,color='red')
        sage: P + Q

    We draw two parametric surfaces and a transparent plane::

        sage: L = plot3d(lambda x,y: 0, (-5,5), (-5,5), color="lightblue", opacity=0.8)
        sage: P = plot3d(lambda x,y: 4 - x^3 - y^2, (-2,2), (-2,2), color='green')
        sage: Q = plot3d(lambda x,y: x^3 + y^2 - 4, (-2,2), (-2,2), color='orange')
        sage: L + P + Q

    We draw the "Sinus" function (water ripple-like surface)::

        sage: x, y = var('x y')
        sage: plot3d(sin(pi*(x^2+y^2))/2,(x,-1,1),(y,-1,1))

    Hill and valley (flat surface with a bump and a dent)::

        sage: x, y = var('x y')
        sage: plot3d( 4*x*exp(-x^2-y^2), (x,-2,2), (y,-2,2))

    An example of a transformation::

        sage: r, phi, z = var('r phi z')
        sage: trans=(r*cos(phi),r*sin(phi),z,z)
        sage: plot3d(cos(r),(r,0,17*pi/2),(phi,0,2*pi),transformation=trans,opacity=0.87).show(aspect_ratio=(1,1,2),frame=False)

    Many more examples of transformations::

        sage: u, v, w = var('u v w')
        sage: rectangular=(u,v,w,w)
        sage: spherical=(w*cos(u)*sin(v),w*sin(u)*sin(v),w*cos(v),w)
        sage: cylindric_radial=(w*cos(u),w*sin(u),v,w)
        sage: cylindric_axial=(v*cos(u),v*sin(u),w,w)
        sage: parabolic_cylindrical=(w*v,(v^2-w^2)/2,u,w)

    Plot a constant function of each of these to get an idea of what it does::

        sage: A = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=rectangular,plot_points=[100,100])
        sage: B = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=spherical,plot_points=[100,100])
        sage: C = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=cylindric_radial,plot_points=[100,100])
        sage: D = plot3d(2,(u,-pi,pi),(v,0,pi),transformation=cylindric_axial,plot_points=[100,100])
        sage: E = plot3d(2,(u,-pi,pi),(v,-pi,pi),transformation=parabolic_cylindrical,plot_points=[100,100])
        sage: @interact
        sage: def _(which_plot=[A,B,C,D,E]):
        ...       show(which_plot)

    Now plot a function::

        sage: g=3+sin(4*u)/2+cos(4*v)/2
        sage: F = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=rectangular,plot_points=[100,100])
        sage: G = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=spherical,plot_points=[100,100])
        sage: H = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=cylindric_radial,plot_points=[100,100])
        sage: I = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=cylindric_axial,plot_points=[100,100])
        sage: J = plot3d(g,(u,-pi,pi),(v,0,pi),transformation=parabolic_cylindrical,plot_points=[100,100])
        sage: @interact
        sage: def _(which_plot=[F, G, H, I, J]):
        ...       show(which_plot)

    TESTS:

    Make sure the transformation plots work::

        sage: show(A + B + C + D + E)
        sage: show(F + G + H + I + J)

    Listing the same plot variable twice gives an error::

        sage: x, y = var('x y')
        sage: plot3d( 4*x*exp(-x^2-y^2), (x,-2,2), (x,-2,2))
        Traceback (most recent call last):
        ...
        ValueError: range variables should be distinct, but there are duplicates
    """
    if transformation is not None:
        from sage.symbolic.callable import is_CallableSymbolicExpression
        # First, determine the parameters for f (from the first item of urange
        # and vrange, preferably).
        if len(urange) == 3 and len(vrange) == 3:
            params = (urange[0], vrange[0])
        elif is_CallableSymbolicExpression(f):
            params = f.variables()
        else:
            # There's no way to find out
            raise ValueError, 'expected 3-tuple for urange and vrange'

        if isinstance(transformation, (tuple, list)):
            transformation = _ArbCoordTrans(transformation[0:3], transformation[3])

        if isinstance(transformation, _CoordTrans):
            R = transformation.to_cartesian(f, params)
            return parametric_plot3d.parametric_plot3d(R, urange, vrange, **kwds)
        else:
            raise ValueError, 'unknown transformation type'
    elif adaptive:
        P = plot3d_adaptive(f, urange, vrange, **kwds)
    else:
        u=fast_float_arg(0)
        v=fast_float_arg(1)
        P=parametric_plot3d.parametric_plot3d((u,v,f), urange, vrange, **kwds)
    P.frame_aspect_ratio([1.0,1.0,0.5])
    return P

def plot3d_adaptive(f, x_range, y_range, color="automatic",
                    grad_f=None,
                    max_bend=.5, max_depth=5, initial_depth=4, num_colors=128, **kwds):
    r"""
    Adaptive 3d plotting of a function of two variables.

    This is used internally by the plot3d command when the option
    ``adaptive=True`` is given.

    INPUT:


    -  ``f`` - a symbolic function or a Python function of
       3 variables.

    -  ``x_range`` - x range of values: 2-tuple (xmin,
       xmax) or 3-tuple (x,xmin,xmax)

    -  ``y_range`` - y range of values: 2-tuple (ymin,
       ymax) or 3-tuple (y,ymin,ymax)

    -  ``grad_f`` - gradient of f as a Python function

    -  ``color`` - "automatic" - a rainbow of num_colors
       colors

    -  ``num_colors`` - (default: 128) number of colors to
       use with default color

    -  ``max_bend`` - (default: 0.5)

    -  ``max_depth`` - (default: 5)

    -  ``initial_depth`` - (default: 4)

    -  ``**kwds`` - standard graphics parameters


    EXAMPLES: We plot `\sin(xy)`::

        sage: from sage.plot.plot3d.plot3d import plot3d_adaptive
        sage: x,y=var('x,y'); plot3d_adaptive(sin(x*y), (x,-pi,pi), (y,-pi,pi), initial_depth=5)
    """
    if initial_depth >= max_depth:
        max_depth = initial_depth

    from sage.plot.misc import setup_for_eval_on_grid
    g, ranges = setup_for_eval_on_grid(f, [x_range,y_range], plot_points=2)
    xmin,xmax = ranges[0][:2]
    ymin,ymax = ranges[1][:2]

    opacity = kwds.get('opacity',1)

    if color == "automatic":
        texture = rainbow(num_colors, 'rgbtuple')
    else:
        if isinstance(color, list):
            texture = color
        else:
            kwds['color'] = color
            texture = Texture(kwds)

    factory = TrivialTriangleFactory()
    plot = TrianglePlot(factory, g, (xmin, xmax), (ymin, ymax), g = grad_f,
                        min_depth=initial_depth, max_depth=max_depth,
                        max_bend=max_bend, num_colors = None)

    P = IndexFaceSet(plot._objects)
    if isinstance(texture, (list, tuple)):
        if len(texture) == 2:
            # do a grid coloring
            xticks = (xmax - xmin)/2**initial_depth
            yticks = (ymax - ymin)/2**initial_depth
            parts = P.partition(lambda x,y,z: (int((x-xmin)/xticks) + int((y-ymin)/yticks)) % 2)
        else:
            # do a topo coloring
            bounds = P.bounding_box()
            min_z = bounds[0][2]
            max_z = bounds[1][2]
            if max_z == min_z:
                span = 0
            else:
                span = (len(texture)-1) / (max_z - min_z)    # max to avoid dividing by 0
            parts = P.partition(lambda x,y,z: int((z-min_z)*span))
        all = []
        for k, G in parts.iteritems():
            G.set_texture(texture[k], opacity=opacity)
            all.append(G)
        P = Graphics3dGroup(all)
    else:
        P.set_texture(texture)

    P.frame_aspect_ratio([1.0,1.0,0.5])
    P._set_extra_kwds(kwds)
    return P

def spherical_plot3d(f, urange, vrange, **kwds):
    """
    Takes a function and plots it in spherical coordinates in the domain specified
    by urange and vrange. This function is equivalent to::

        sage: var('r,u,u')
        sage: T = (r*cos(u)*sin(v), r*sin(u)*sin(v), r*cos(v), r)
        sage: plot3d(f, urange, vrange, transformation=T)

    INPUT:

    - ``f`` - a symbolic expression or function of two variables.

    - ``urange`` - a 3-tuple (u, u_min, u_max), the domain of the azimuth variable.

    - ``vrange`` - a 3-tuple (v, v_min, v_max), the domain of the inclination variable.

    EXAMPLES:

    A sphere of radius 2::

        sage: spherical_plot3d(2,(x,0,2*pi),(y,0,pi))

    A drop of water::

        sage: spherical_plot3d(e^-y,(x,0,2*pi),(y,0,pi),opacity=0.5).show(frame=False)

    An object similar to a heart::

        sage: spherical_plot3d((2+cos(2*x))*(y+1),(x,0,2*pi),(y,0,pi),rgbcolor=(1,.1,.1))

    Some random figures:

    ::

        sage: spherical_plot3d(1+sin(5*x)/5,(x,0,2*pi),(y,0,pi),rgbcolor=(1,0.5,0),plot_points=(80,80),opacity=0.7)

    ::

        sage: spherical_plot3d(1+2*cos(2*y),(x,0,3*pi/2),(y,0,pi)).show(aspect_ratio=(1,1,1))
    """
    return plot3d(f, urange, vrange, transformation=Spherical('r', ['phi', 'theta']), **kwds)

def cylindrical_plot3d(f, urange, vrange, **kwds):
    """
    Takes a function and plots it in cylindrical coordinates in the domain specified
    by urange and vrange. This command is equivalent to::

        sage: var('r,u,v')
        sage: T = (r*cos(u), r*sin(u), v, r)
        sage: plot3d(f, urange, vrange, transformation=T)

    INPUT:

    - ``f`` - a symbolic expression or function of two variables.

    - ``urange`` - a 3-tuple (u, u_min, u_max), the domain of the azimuth variable.

    - ``vrange`` - a 3-tuple (v, v_min, v_max), the domain of the elevation (z) variable.

    EXAMPLES:

    A portion of a cylinder of radius 2::

        sage: var('fi,z')
        sage: cylindrical_plot3d(2,(fi,0,3*pi/2),(z,-2,2))

    Some random figures:

    ::

        sage: cylindrical_plot3d(cosh(z),(fi,0,2*pi),(z,-2,2))

    ::

        sage: cylindrical_plot3d(e^(-z^2)*(cos(4*fi)+2)+1,(fi,0,2*pi),(z,-2,2),plot_points=[80,80]).show(aspect_ratio=(1,1,1))
    """
    return plot3d(f, urange, vrange, transformation=Cylindrical('rho', ['phi', 'z']), **kwds)

def axes(scale=1, radius=None, **kwds):
    if radius is None:
        radius = scale/100.0
    return Graphics3dGroup([arrow3d((0,0,0),(scale,0,0), radius, **kwds),
                            arrow3d((0,0,0),(0,scale,0), radius, **kwds),
                            arrow3d((0,0,0),(0,0,scale), radius, **kwds)])
