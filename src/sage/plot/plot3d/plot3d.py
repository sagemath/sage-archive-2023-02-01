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

class TrivialTriangleFactory:
    def triangle(self, a, b, c, color = None):
        return [a,b,c]
    def smooth_triangle(self, a, b, c, da, db, dc, color = None):
        return [a,b,c]

import parametric_plot3d
def plot3d(f, urange, vrange, adaptive=False, **kwds):
    """
    INPUT:


    -  ``f`` - a symbolic expression or function of 2
       variables

    -  ``urange`` - a 2-tuple (u_min, u_max) or a 3-tuple
       (u, u_min, u_max)

    -  ``vrange`` - a 2-tuple (v_min, v_max) or a 3-tuple
       (v, v_min, v_max)

    -  ``adaptive`` - (default: False) whether to use
       adaptive refinement to draw the plot (slower, but may look better)

    -  ``mesh`` - bool (default: False) whether to display
       mesh grid lines

    -  ``dots`` - bool (default: False) whether to display
       dots at mesh grid points

    -  ``plot_points`` - (default: "automatic") initial number of sample
       points in each direction; an integer or a pair of integers

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

    TESTS: Listing the same plot variable twice gives an error.

    ::

        sage: x, y = var('x y')
        sage: plot3d( 4*x*exp(-x^2-y^2), (x,-2,2), (x,-2,2))
        Traceback (most recent call last):
        ...
        ValueError: range variables should be distinct, but there are duplicates
    """
    if adaptive:
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


def axes(scale=1, radius=None, **kwds):
    if radius is None:
        radius = scale/100.0
    return Graphics3dGroup([arrow3d((0,0,0),(scale,0,0), radius, **kwds),
                            arrow3d((0,0,0),(0,scale,0), radius, **kwds),
                            arrow3d((0,0,0),(0,0,scale), radius, **kwds)])
