r"""
Functions for 3d plotting with the new graphcis model.

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
    -- Robert Bradshaw 2007-08: inital version of this file

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

class TrivialTriangleFactory:
    def triangle(self, a, b, c, color = None):
        return [a,b,c]
    def smooth_triangle(self, a, b, c, da, db, dc, color = None):
        return [a,b,c]


def plot3d(f,(xmin,xmax),(ymin,ymax),texture,grad_f=None,
              max_bend=.7,max_depth=5,initial_depth=3, num_colors=None):

    factory = TrivialTriangleFactory()
    plot = TrianglePlot(factory, f, (xmin, xmax), (ymin, ymax), g = grad_f,
                         min_depth=initial_depth, max_depth=max_depth, max_bend=max_bend, num_colors = num_colors)

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
            span = len(texture) / (max_z - min_z)
            parts = P.partition(lambda x,y,z: int((z-min_z)*span))
        all = []
        for k, G in parts.iteritems():
            G.set_texture(texture[k])
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