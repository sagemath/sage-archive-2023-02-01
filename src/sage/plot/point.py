"""

TESTS:
    sage: E = EllipticCurve('37a')
    sage: P = E(0,0)
    sage: def get_points(n): return sum([point(list(i*P)[:2], pointsize=3) for i in range(-n,n) if i != 0 and (i*P)[0] < 3])
    sage: sum([get_points(15*n).plot3d(z=n) for n in range(1,10)])
"""

from sage.plot.primitive import GraphicPrimitive_xydata
#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
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
from sage.plot.misc import options, rename_keyword, to_mpl_color

class Point(GraphicPrimitive_xydata):
    """
    Primitive class that initializes the
    point graphics type

    """
    def __init__(self, xdata, ydata, options):
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive_xydata.__init__(self, options)

    def _allowed_options(self):
        return {'alpha':'How transparent the line is.',
                'pointsize': 'How big the point is.',
                'faceted': 'If True color the edge of the point.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

    def _plot3d_options(self, options=None):
        if options == None:
            options = dict(self.options())
        options_3d = {}
        if 'pointsize' in options:
            options_3d['size'] = options['pointsize']
            del options['pointsize']
        if 'faceted' in options:
            if options['faceted']:
                raise NotImplementedError, "No 3d faceted points."
            del options['faceted']
        options_3d.update(GraphicPrimitive_xydata._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: P = E(0,0)
            sage: def get_points(n):
            ...     return sum([point((i*P)[0:2], pointsize=3) for i in range(-n,n) if i != 0 and (i*P)[0] < 3])
            sage: sum([get_points(15*n).plot3d(z=n) for n in range(1,10)])
        """
        from sage.plot.plot3d.base import Graphics3dGroup
        from sage.plot.plot3d.shapes2 import point3d
        options = self._plot3d_options()
        options.update(kwds)
        all = [point3d([(x, y, 0) for x, y in zip(self.xdata, self.ydata)], **options)]
        if len(all) == 1:
            return all[0]
        else:
            return Graphics3dGroup(all)

    def _repr_(self):
        return "Point set defined by %s point(s)"%len(self.xdata)

    def __getitem__(self, i):
        return self.xdata[i], self.ydata[i]

    def _render_on_subplot(self,subplot):
        r"""
        TESTS:
        We check to make sure that \#2076 is fixed by verifying all
        the points are red.
            sage: point(((1,1), (2,2), (3,3)), rgbcolor=hue(1), pointsize=30)
        """
        options = self.options()

        #Convert the color to a hex string so that the scatter
        #method does not interpret it as a list of 3 floating
        #point color specifications when there are
        #three points. This is mentioned in the matplotlib 0.98
        #documentation and fixes \#2076
        from matplotlib.colors import rgb2hex
        c = rgb2hex(to_mpl_color(options['rgbcolor']))

        a = float(options['alpha'])
        s = int(options['pointsize'])
        faceted = options['faceted'] #faceted=True colors the edge of point
        scatteroptions={}
        if not faceted: scatteroptions['edgecolors'] = 'none'

        subplot.scatter(self.xdata, self.ydata, s=s, c=c, alpha=a, **scatteroptions)

def point(points, **kwds):
    """
    Returns either a 2-dimensional or 3-dimensional point or sum of points.

    INPUT:
        points -- either a single point (as a tuple) or a list of points.

    For information regarding additional arguments, see either point2d?
    or point3d?.

    EXAMPLES:
        sage: point((1,2))
        sage: point((1,2,3))
        sage: point([(0,0), (1,1)])
        sage: point([(0,0,1), (1,1,1)])
    """
    try:
        return point2d(points, **kwds)
    except (ValueError, TypeError):
        from sage.plot.plot3d.shapes2 import point3d
        return point3d(points, **kwds)

@options(alpha=1, pointsize=10, faceted=False, rgbcolor=(0,0,1))
def point2d(points, **options):
    r"""
    A point of size `pointsize' defined by point = $(x,y)$.
    Point takes either a single tuple of coordinates or a list of tuples.

    Type \code{point2d.options} to see all options.

    EXAMPLES:
        A purple point from a single tuple or coordinates:
        sage: point((0.5, 0.5), rgbcolor=hue(0.75))

        Here are some random larger red points, given as a list of tuples
        sage: point(((0.5, 0.5), (1, 2), (0.5, 0.9), (-1, -1)), rgbcolor=hue(1), pointsize=30)

    """
    from sage.plot.plot import xydata_from_point_list, Graphics
    xdata, ydata = xydata_from_point_list(points)
    g = Graphics()
    g.add_primitive(Point(xdata, ydata, options))
    return g

points = point
