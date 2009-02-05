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
from primitive import GraphicPrimitive
from sage.plot.misc import options, rename_keyword, to_mpl_color
from math import sin, cos, pi

class Circle(GraphicPrimitive):
    """
    Circle graphics primitive.
    """
    def __init__(self, x, y, r, options):
        self.x = float(x)
        self.y = float(y)
        self.r = float(r)
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES:
            sage: p = circle((3, 3), 1)
            sage: d = p.get_minmax_data()
            sage: d['xmin']
            2.0
            sage: d['ymin']
            2.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data([self.x - self.r, self.x + self.r],
                           [self.y - self.r, self.y + self.r],
                           dict=True)

    def _allowed_options(self):
        return {'alpha':'How transparent the line is.',
                'fill': 'Whether or not to fill the polygon.',
                'thickness':'How thick the border of the polygon is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'zorder':'The layer level in which to draw',
                'linestyle':"The style of the line, which is one of 'dashed', 'dotted', 'solid', 'dashdot'."}

    def _repr_(self):
        return "Circle defined by (%s,%s) with r=%s"%(self.x, self.y, self.r)

    def _render_on_subplot(self, subplot):
        import matplotlib.patches as patches
        options = self.options()
        p = patches.Circle((float(self.x), float(self.y)), float(self.r))
        p.set_linewidth(float(options['thickness']))
        p.set_fill(options['fill'])
        a = float(options['alpha'])
        p.set_alpha(a)
        c = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(c)
        p.set_facecolor(c)
        p.set_linestyle(options['linestyle'])
        subplot.add_patch(p)

    def plot3d(self, **kwds):
        """
        EXAMPLES:
            sage: circle((0,0), 1).plot3d()
            sage: sum([circle((random(),random()), random()).plot3d(z=random()) for _ in range(20)])
        """
        options = dict(self.options())
        fill = options['fill']
        del options['fill']
        del options['linestyle']
        n = 50
        dt = float(2*pi/n)
        x, y, r = self.x, self.y, self.r
        xdata = [x+r*cos(t*dt) for t in range(n+1)]
        ydata = [y+r*sin(t*dt) for t in range(n+1)]
        if fill:
            from polygon import Polygon
            return Polygon(xdata, ydata, options).plot3d()
        else:
            from line import Line
            return Line(xdata, ydata, options).plot3d()


@options(alpha=1, fill=False, thickness=1, rgbcolor=(0,0,1), linestyle='solid')
def circle(center, radius, **options):
    """
    Return a circle at a point = $(x,y)$ with radius = $r$.
    Type \code{circle.options} to see all options.

    circle(center, radius, **kwds)

    INPUT:
        center -- a 2-tuple (x,y)
        radius -- a positive number
        alpha -- default: 1
        fill -- default: False
        thickness -- default: 1
        rgbcolor -- default: (0,0,0)
        linestyle -- default: 'solid' (2d plotting only)

    EXAMPLES:
        sage: c = circle((1,1), 1, rgbcolor=(1,0,0))
        sage: c

    To correct the apect ratio of certain graphics, it is necessary
    to show with a `\code{figsize}' of square dimensions.

        sage: c.show(figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)

    Here we make an more complicated plot with many circles of different colors

        sage: g = Graphics()
        sage: step=6; ocur=1/5; paths=16;
        sage: PI = math.pi    # numerical for speed -- fine for graphics
        sage: for r in range(1,paths+1):
        ...       for x,y in [((r+ocur)*math.cos(n), (r+ocur)*math.sin(n)) for n in srange(0, 2*PI+PI/step, PI/step)]:
        ...           g += circle((x,y), ocur, rgbcolor=hue(r/paths))
        ...       rnext = (r+1)^2
        ...       ocur = (rnext-r)-ocur
        ...
        sage: g.show(xmin=-(paths+1)^2, xmax=(paths+1)^2, ymin=-(paths+1)^2, ymax=(paths+1)^2, figsize=[6,6])

    """
    from sage.plot.plot import Graphics
    g = Graphics()
    g.add_primitive(Circle(center[0], center[1], radius, options))
    return g
