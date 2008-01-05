r"""
Functions for 3d plotting with the new graphics model.

EXAMPLES:
    sage: from sage.plot.plot3d.plot3d import plot3d, axes
    sage: from sage.plot.plot import rainbow

    sage: def f(x,y):
    ...       return math.sin(y*y+x*x)/math.sqrt(x*x+y*y+.0001)
    ...
    sage: P = plot3d(f,(-3,3),(-3,3), rainbow(60, 'rgbtuple'), max_bend=.1, max_depth=15)
    sage: P.show()

    sage: def f(x,y):
    ...       return math.exp(x/5)*math.sin(y)
    ...
    sage: P = plot3d(f,(-5,5),(-5,5), ['red','yellow'])
    sage: S = P + axes(6, color='black')
    sage: S.show()

AUTHOR:
    -- Tom Boothby: adaptive refinement triangles
    -- Joshua Kantor: adaptive refinement triangles
    -- Robert Bradshaw 2007-08: initial version of this file

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
from shapes import Arrow
from base import Graphics3dGroup
from sage.plot.plot import rainbow

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
        sage: show(plot3d(lambda x, y: x^2 + y^2, (-2,2), (-2,2)))

    We plot some 3d symbolic functions:
        sage: var('x,y')
        sage: show(plot3d(x^2 + y^2, (x,-2,2), (y,-2,2)))
        sage: show(plot3d(sin(x*y), (x, -pi, pi), (y, -pi, pi)))

    We draw two parametric surfaces and a transparent plane:
        sage: L = plot3d(lambda x,y: 0, (-5,5), (-5,5), texture=Texture("lightblue", opacity=0.8))
        sage: P = plot3d(lambda x,y: 4 - x^3 - y^2, (-2,2), (-2,2), texture=Texture('green'))
        sage: Q = plot3d(lambda x,y: x^3 + y^2 - 4, (-2,2), (-2,2), texture=Texture('orange'))
        sage: show(L + P + Q)
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

def plot3d_rainbow(f,(xmin,xmax),(ymin,ymax), texture=None, opacity=1, grad_f=None,
                   max_bend=.5, max_depth=5, initial_depth=4, num_colors=None):
    """
    EXAMPLES:
    """
    if initial_depth >= max_depth:
        max_depth = initial_depth
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)

    # Check if f has a fast float evaluation
    #try:
    #    f = f.fast_float_function()
    #except AttributeError:
    #    # Nope -- no prob.
    #    pass
    if texture is None:
        texture = rainbow(128, 'rgbtuple')

    factory = TrivialTriangleFactory()
    plot = TrianglePlot(factory, f, (xmin, xmax), (ymin, ymax), g = grad_f,
                        min_depth=initial_depth, max_depth=max_depth,
                        max_bend=max_bend, num_colors = num_colors)

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
        return Graphics3dGroup(all)
    else:
        P.set_texture(texture)
        return P

def axes(scale=1, radius=None, **kwds):
    if radius is None:
        radius = scale/100.0
    return Graphics3dGroup([Arrow((0,0,0),(scale,0,0), radius, **kwds),
                            Arrow((0,0,0),(0,scale,0), radius, **kwds),
                            Arrow((0,0,0),(0,0,scale), radius, **kwds)])
