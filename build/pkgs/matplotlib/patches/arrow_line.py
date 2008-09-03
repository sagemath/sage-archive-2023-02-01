r"""
A matplotlib subclass to draw an arrowhead on a line.

AUTHORS:
    -- Jason Grout (2008-08-19): initial version
"""

############################################################################
# Copyright (c) 2008, Jason Grout <jason-sage@creativetrax.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
############################################################################

import matplotlib
from matplotlib.path import Path
from matplotlib.lines import Line2D
import math
import matplotlib.cbook
from matplotlib.transforms import Affine2D

class ArrowLine(Line2D):
    """
    A matplotlib subclass to draw an arrowhead on a line.

    INPUT:
        arrow -- (default '>') What sort of arrow to draw.  Currently the options are just '>', a triangular arrow.
        arrowsize -- (default 8)
        arrowedge -- (default False) whether or not to draw an edge around the arrowhead.  If an edge is drawn, the arrowhead does not precisely point to the endpoint (the edge goes past the endpoint).
        arrowheadlength -- (default arrowsize)
        arrowheadwidth -- (default arrowsize)
        arrowfacecolor, arrowedgecolor -- the colors for the arrowhead face and edge.
        arrowedgewidth -- (default 4) the width of the edge around the arrowhead
        arrowshorten -- (default 0) the number of points by which to shorten the arrow.  The shortened arrow is centered on the old arrow.  This allows an arrow, for example to easily go between the boundaries of a circle by specifying the endpoints as the center of the circles and the arrowshorten as twice the radii of the circles.

    EXAMPLE:
        sage: import pylab
        sage: fig = pylab.figure()
        sage: ax = fig.add_subplot(111, autoscale_on=False)
        sage: t = [-1,2]
        sage: s = [0,-1]
        sage: line = ArrowLine(t, s, color='b', ls='-', lw=2, arrow='>', arrowsize=20)
        sage: ax.add_line(line)
        sage: ax.set_xlim(-3,3)
        (-3, 3)
        sage: ax.set_ylim(-3,3)
        (-3, 3)
        sage: pylab.savefig('test.png')

    """

    arrows = {'>' : '_draw_triangle_arrow'}

    def __init__(self, *args, **kwargs):
        """Initialize the line and arrow."""
        self._arrow = kwargs.pop('arrow', None)
        self._arrowsize = kwargs.pop('arrowsize', 8)
        self._arrowedge = kwargs.pop('arrowedge',False)
        self._arrowedgecolor = kwargs.pop('arrowedgecolor', kwargs.get('rgbcolor', 'k'))
        self._arrowfacecolor = kwargs.pop('arrowfacecolor', kwargs.get('rgbcolor', 'k'))
        self._arrowedgewidth = kwargs.pop('arrowedgewidth', 4)
        self._arrowheadwidth = kwargs.pop('arrowheadwidth', self._arrowsize)
        self._arrowheadlength = kwargs.pop('arrowheadlength', self._arrowsize)
        self._arrowshorten = kwargs.pop('arrowshorten', 0)
        Line2D.__init__(self, *args, **kwargs)



    def draw(self, renderer):
        """Draw the line and arrowhead using the passed renderer."""
        if self._invalid:
            self.recache()
        renderer.open_group('arrowline2d')
        if not self._visible: return

        orig_t = self.get_transform()
        points = self.get_xydata()
        pixel_vector = (orig_t.transform_point(points[1]) - orig_t.transform_point(points[0]))
        pixel_length=math.sqrt(sum(pixel_vector**2))
        clip_fraction = renderer.points_to_pixels(self._arrowshorten)/pixel_length
        head_clip_fraction = renderer.points_to_pixels(self._arrowshorten+self._arrowheadlength*0.8)/pixel_length

        real_vector = (points[1]-points[0])
        translate = (real_vector)*(clip_fraction)*0.5

        to_origin = Affine2D().translate(-points[0][0],-points[0][1])
        # from_origin also includes the translation which moves the
        # arrow so that it is centered on its former self.
        from_origin = Affine2D().translate(translate[0],translate[1]) \
            + Affine2D().translate(points[0][0],points[0][1])

        # Now pull the line back so it doesn't stick through the arrow:
        self.set_transform(to_origin \
                               + Affine2D().scale(1.0-head_clip_fraction) \
                               + from_origin \
                               + orig_t)
        self.recache()

        Line2D.draw(self, renderer)

        self.set_transform(to_origin \
                               + Affine2D().scale(1.0-clip_fraction) \
                               + from_origin \
                               + orig_t)

        self.recache()

        if self._arrow is not None:
            gc = renderer.new_gc()
            self._set_gc_clip(gc)
            if self._arrowedge:
                gc.set_foreground(self._arrowedgecolor)
                gc.set_linewidth(self._arrowedgewidth)
            else:
                gc.set_linewidth(0) # set linewidth to zero to get a true triangle ending at the endpoint.
            gc.set_alpha(self._alpha)
            funcname = self.arrows.get(self._arrow, '_draw_nothing')
            if funcname != '_draw_nothing':
                tpath, affine = self._transformed_path.get_transformed_points_and_affine()
                arrowFunc = getattr(self, funcname)
                arrowFunc(renderer, gc, tpath, affine.frozen())

        renderer.close_group('arrowline2d')

        self.set_transform(orig_t)



    _arrow_path = Path([[0.0, 0.0], [-1.0, 1.0], [-1.0, -1.0], [0.0, 0.0]], codes=[Path.MOVETO, Path.LINETO,Path.LINETO, Path.CLOSEPOLY])
    def _draw_triangle_arrow(self, renderer, gc, path, path_trans):
        """Draw a triangular arrow."""
        segment = [i[0] for i in path.iter_segments()][-2:]
        startx,starty = path_trans.transform_point(segment[0])
        endx,endy = path_trans.transform_point(segment[1])
        angle = math.atan2(endy-starty, endx-startx)
        halfwidth = 0.5*renderer.points_to_pixels(self._arrowheadwidth)
        length = renderer.points_to_pixels(self._arrowheadlength)
        transform = Affine2D().scale(length,halfwidth).rotate(angle).translate(endx,endy)

        rgbFace = self._get_rgb_arrowface()
        renderer.draw_path(gc, self._arrow_path, transform, rgbFace)

    def _get_rgb_arrowface(self):
        facecolor = self._arrowfacecolor
        if matplotlib.cbook.is_string_like(facecolor) and facecolor.lower()=='none':
            rgbFace = None
        else:
            rgbFace = matplotlib.colors.colorConverter.to_rgb(facecolor)
        return rgbFace



