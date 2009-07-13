"""
Circles
"""
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
from sage.plot.misc import options, rename_keyword
from sage.plot.colors import to_mpl_color
from math import sin, cos, pi

class Circle(GraphicPrimitive):
    """
    Primitive class for the Circle graphics type.  See circle? for information
    about actually plotting circles.

    INPUT:

    - x - `x`-coordinate of center of Circle

    - y - `y`-coordinate of center of Circle

    - r - radius of Circle object

    - options - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via ``circle``::

        sage: from sage.plot.circle import Circle
        sage: C = Circle(2,3,5,{'zorder':2})
        sage: C
        Circle defined by (2.0,3.0) with r=5.0
        sage: C.options()['zorder']
        2
        sage: C.r
        5.0

    TESTS:

    We test creating a circle::

        sage: C = circle((2,3), 5)
    """
    def __init__(self, x, y, r, options):
        """
        Initializes base class Circle.

        EXAMPLES::

            sage: C = circle((2,3), 5, edgecolor='red', alpha=.5, fill=True)
            sage: C[0].x
            2.0
            sage: C[0].r
            5.0
            sage: C[0].options()['edgecolor']
            'red'
            sage: C[0].options()['alpha']
            0.500000000000000
        """
        self.x = float(x)
        self.y = float(y)
        self.r = float(r)
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

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
        """
        Return the allowed options for the Circle class.

        EXAMPLES::

            sage: p = circle((3, 3), 1)
            sage: p[0]._allowed_options()['alpha']
            'How transparent the figure is.'
            sage: p[0]._allowed_options()['facecolor']
            '2D only: The color of the face as an rgb tuple.'
        """
        return {'alpha':'How transparent the figure is.',
                'fill': 'Whether or not to fill the circle.',
                'thickness':'How thick the border of the circle is.',
                'edgecolor':'2D only: The color of the edge as an rgb tuple.',
                'facecolor':'2D only: The color of the face as an rgb tuple.',
                'rgbcolor':'The color (edge and face) as an rgb tuple.',
                'hue':'The color given as a hue.',
                'zorder':'2D only: The layer level in which to draw',
                'linestyle':"2D only: The style of the line, which is one of 'dashed', 'dotted', 'solid', 'dashdot'."}

    def _repr_(self):
        """
        String representation of Circle primitive.

        EXAMPLES::

            sage: C = circle((2,3), 5)
            sage: c = C[0]; c
            Circle defined by (2.0,3.0) with r=5.0
        """
        return "Circle defined by (%s,%s) with r=%s"%(self.x, self.y, self.r)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: C = circle((2,pi), 2, edgecolor='black', facecolor='green', fill=True)
        """
        import matplotlib.patches as patches
        options = self.options()
        p = patches.Circle((float(self.x), float(self.y)), float(self.r))
        p.set_linewidth(float(options['thickness']))
        p.set_fill(options['fill'])
        a = float(options['alpha'])
        p.set_alpha(a)
        ec = to_mpl_color(options['edgecolor'])
        fc = to_mpl_color(options['facecolor'])
        if 'rgbcolor' in options:
            ec = fc = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(ec)
        p.set_facecolor(fc)
        p.set_linestyle(options['linestyle'])
        subplot.add_patch(p)

    def plot3d(self, z=0, **kwds):
        """
        Plots a 2D circle (actually a 50-gon) in 3D,
        with default height zero.

        INPUT:


        -  ``z`` - optional 3D height above `xy`-plane.

        EXAMPLES:

            sage: circle((0,0), 1).plot3d()

        This example uses this method implicitly, but does not pass
        the optional parameter z to this method::

            sage: sum([circle((random(),random()), random()).plot3d(z=random()) for _ in range(20)])

        These examples are explicit, and pass z to this method::

            sage: C = circle((2,pi), 2, hue=.8, alpha=.3, fill=True)
            sage: c = C[0]
            sage: d = c.plot3d(z=2)
            sage: d.texture.opacity
            0.300000000000000

        ::

            sage: C = circle((2,pi), 2, hue=.8, alpha=.3)
            sage: c = C[0]
            sage: d = c.plot3d(z=2)
            sage: d.jmol_repr(d.testing_render_params())[0][-1]
            'color $line_1 translucent 0.7 [204,0,255]'
        """
        options = dict(self.options())
        fill = options['fill']
        del options['fill']
        del options['linestyle']
        del options['edgecolor']
        del options['facecolor']
        del options['zorder']
        n = 50
        dt = float(2*pi/n)
        x, y, r = self.x, self.y, self.r
        xdata = [x+r*cos(t*dt) for t in range(n+1)]
        ydata = [y+r*sin(t*dt) for t in range(n+1)]
        if fill:
            from polygon import Polygon
            return Polygon(xdata, ydata, options).plot3d(z)
        else:
            from line import Line
            return Line(xdata, ydata, options).plot3d().translate((0,0,z))


@options(alpha=1, fill=False, thickness=1, edgecolor='black', facecolor='red', linestyle='solid', zorder=5)
def circle(center, radius, **options):
    """
    Return a circle at a point center = `(x,y)` (or `(x,y,z)` and
    parallel to the `xy`-plane) with radius = `r`.  Type
    ``circle.options`` to see all options.

    OPTIONS:

    - ``alpha`` - default: 1

    - ``fill`` - default: False

    - ``thickness`` - default: 1

    - ``rgbcolor`` - default: (0,0,0)

    - ``linestyle`` - default: 'solid' (2D plotting only)

    - ``edgecolor`` - default: 'black' (2D plotting only)

    - ``facecolor`` - default: 'red' (2D plotting only, useful only
      if Fill=True)

    EXAMPLES::

        sage: c = circle((1,1), 1, rgbcolor=(1,0,0))
        sage: c

    We can also use this command to plot three-dimensional circles parallel
    to the `xy`-plane.

        sage: c = circle((1,1,3), 1, rgbcolor=(1,0,0))
        sage: c
        sage: type(c)
        <class 'sage.plot.plot3d.base.TransformGroup'>

    To correct the aspect ratio of certain graphics, it is necessary
    to show with a ``figsize`` of square dimensions::

        sage: c.show(figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)

    Here we make a more complicated plot, with many circles of different colors::

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

    Note that the ``rgbcolor`` option overrides the other coloring options.
    This produces red fill in a blue circle::

        sage: C = circle((2,3), 1, fill=True, edgecolor='blue')

    This produces an all-green filled circle::

        sage: C = circle((2,3), 1, fill=True, edgecolor='blue', rgbcolor='green')

    The option ``hue`` overrides *all* other options, so be careful with its use.
    This produces a purplish filled circle::

        sage: C = circle((2,3), 1, fill=True, edgecolor='blue', rgbcolor='green', hue=.8)

    TESTS:

    We cannot currently plot circles in more than three dimensions::

        sage: c = circle((1,1,1,1), 1, rgbcolor=(1,0,0))
        Traceback (most recent call last):
        ...
        ValueError: The center of a plotted circle should have two or three coordinates.
    """
    from sage.plot.plot import Graphics
    g = Graphics()
    g.add_primitive(Circle(center[0], center[1], radius, options))
    if len(center)==2:
        return g
    elif len(center)==3:
        return g[0].plot3d(z=center[2])
    else:
        raise ValueError, 'The center of a plotted circle should have two or three coordinates.'
