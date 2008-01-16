r"""
Functions for 3d plotting with the new graphics model.

EXAMPLES:
    sage: def f(x,y):
    ...       return math.sin(y*y+x*x)/math.sqrt(x*x+y*y+.0001)
    ...
    sage: P = plot3d_adaptive(f,(-3,3),(-3,3), color=rainbow(60, 'rgbtuple'), max_bend=.1, max_depth=15)
    sage: P.show()

    sage: def f(x,y):
    ...       return math.exp(x/5)*math.sin(y)
    ...
    sage: P = plot3d_adaptive(f,(-5,5),(-5,5), color=['red','yellow'])
    sage: from sage.plot.plot3d.plot3d import axes
    sage: S = P + axes(6, color='black')
    sage: S.show()

We plot "cape man":
    sage: S = sphere(size=.5, color='yellow')

    sage: from sage.plot.plot3d.shapes import Cone
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)

    sage: S += sphere((.45,-.1,.15), size=.1, color='white') + sphere((.51,-.1,.17), size=.05, color='black')
    sage: S += sphere((.45, .1,.15),size=.1, color='white') + sphere((.51, .1,.17), size=.05, color='black')
    sage: S += sphere((.5,0,-.2),size=.1, color='yellow')
    sage: def f(x,y): return math.exp(x/5)*math.cos(y)
    sage: P = plot3d_adaptive(f,(-5,5),(-5,5), ['red','yellow'], max_depth=10)
    sage: cape_man = P.scale(.2) + S.translate(1,0,0)
    sage: cape_man

AUTHOR:
    -- Tom Boothby: adaptive refinement triangles
    -- Joshua Kantor: adaptive refinement triangles
    -- Robert Bradshaw 2007-08: initial version of this file
    -- William Stein 2007-12, 2008-01: improving 3d plotting

TODO:
    -- smooth triangles
"""

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


from sage.plot.tri_plot import TrianglePlot
from index_face_set import IndexFaceSet
from shapes import arrow3d
from base import Graphics3dGroup
from sage.plot.plot import rainbow
from texture import Texture, is_Texture

class TrivialTriangleFactory:
    def triangle(self, a, b, c, color = None):
        return [a,b,c]
    def smooth_triangle(self, a, b, c, da, db, dc, color = None):
        return [a,b,c]

import parametric_plot3d
def plot3d(f, urange, vrange, **kwds):
    """
    EXAMPLES:
    We plot a 3d function defined as a Python function:
        sage: plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2))

    We plot some 3d symbolic functions:
        sage: x, y = var('x,y')
        sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
        sage: plot3d(sin(x*y), (x, -pi, pi), (y, -pi, pi))

    Two wobby translucent planes:
        sage: x,y = var('x,y')
        sage: P = plot3d(x+y+sin(x*y), (x,-10,10),(y,-10,10), opacity=0.87)
        sage: Q = plot3d(x-2*y-cos(x*y),(x,-10,10),(y,-10,10),opacity=0.3,color='red')
        sage: P + Q

    We draw two parametric surfaces and a transparent plane:
        sage: L = plot3d(lambda x,y: 0, (-5,5), (-5,5), color="lightblue", opacity=0.8)
        sage: P = plot3d(lambda x,y: 4 - x^3 - y^2, (-2,2), (-2,2), color='green')
        sage: Q = plot3d(lambda x,y: x^3 + y^2 - 4, (-2,2), (-2,2), color='orange')
        sage: L + P + Q
    """
    if len(urange) == 2:
        w = (lambda u,v: u, lambda u,v: v, f)
    else:
        u = urange[0]
        v = vrange[0]
        w = (u, v, f)
    P = parametric_plot3d.parametric_plot3d(w, urange, vrange, **kwds)
    P.frame_aspect_ratio([1.0,1.0,0.5])
    return P

def plot3d_adaptive(f, x_range, y_range, color="automatic",
                    grad_f=None,
                    max_bend=.5, max_depth=5, initial_depth=4, num_colors=128, **kwds):
    """
    Adaptive 3d plotting of a function of two variables.

    INPUTS:
        f -- a symbolic function or a Python function of 3 variables.
        x_range -- x range of values: 2-tuple (xmin, xmax) or 3-tuple (x,xmin,xmax)
        y_range -- y range of values: 2-tuple (ymin, ymax) or 3-tuple (y,ymin,ymax)
        grad_f -- gradient of f as a Python function
        color -- "automatic" -- a rainbow of num_colors colors
        num_colors -- (default: 128) number of colors to use with default color
        max_bend -- (default: 0.5)
        max_depth -- (default: 5)
        initial_depth -- (default: 4)
        **kwds -- standard graphics parameters

    EXAMPLES:
    We plot sin(xy):
        sage: x,y=var('x,y'); plot3d_adaptive(sin(x*y), (x,-pi,pi), (y,-pi,pi), initial_depth=5)

    """
    if not (isinstance(x_range, (tuple, list)) and isinstance(y_range,(tuple,list)) and len(x_range) == len(y_range) and len(x_range) in [2,3]):
        raise TypeError, "x_range and y_range must both be a 2 or 3-tuple"
    if len(x_range) == 3:
        # symbolic case
        x, xmin, xmax = x_range
        y, ymin, ymax = y_range
        def g(xx,yy):
            return float(f.subs({x:xx, y:yy}))

    else:
        xmin, xmax = x_range
        ymin, ymax = y_range

        g = f

    if initial_depth >= max_depth:
        max_depth = initial_depth
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)

    # Check if g has a fast float evaluation
    #try:
    #    g = g.fast_float_function()
    #except AttributeError:
    #    # Nope -- no prob.
    #    pass
    if kwds.has_key('opacity'):
        opacity = kwds['opacity']
    else:
        opacity = 1
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
                span = len(texture) / (max_z - min_z)    # max to avoid dividing by 0
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
