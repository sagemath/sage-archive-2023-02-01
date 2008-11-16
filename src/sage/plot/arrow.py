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
from sage.plot.primitive import GraphicPrimitive
from sage.plot.misc import options, rename_keyword, to_mpl_color

class Arrow(GraphicPrimitive):
    """
    Primitive class that initializes the arrow graphics type

    EXAMPLES:
    We create an arrow graphics object, then take the 0th entry
    in it to get the actual Arrow graphics primitive:
        sage: P = arrow((0,1), (2,3))[0]
        sage: type(P)
        <class 'sage.plot.arrow.Arrow'>
        sage: P
        Arrow from (0.0,1.0) to (2.0,3.0)
    """
    def __init__(self, xtail, ytail, xhead, yhead, options):
        """
        Create an arrow graphics primitive.

        EXAMPLES:
            sage: from sage.plot.arrow import Arrow
            sage: Arrow(0,0,2,3,{})
            Arrow from (0.0,0.0) to (2.0,3.0)
        """
        self.xtail = float(xtail)
        self.xhead = float(xhead)
        self.ytail = float(ytail)
        self.yhead = float(yhead)
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a bounding box for this arrow.

        EXAMPLES:
            sage: d = arrow((1,1), (5,5)).get_minmax_data()
            sage: d['xmin']
            1.0
            sage: d['xmax']
            5.0
        """
        return {'xmin': min(self.xtail, self.xhead),
                'xmax': max(self.xtail, self.xhead),
                'ymin': min(self.ytail, self.yhead),
                'ymax': max(self.ytail, self.yhead)}


    def _allowed_options(self):
        """
        Return the dictionary of allowed options for the arrow graphics primitive.

        EXAMPLES:
             sage: from sage.plot.arrow import Arrow
             sage: list(sorted(Arrow(0,0,2,3,{})._allowed_options().iteritems()))
             [('arrowshorten', 'The length in points to shorten the arrow.'),
             ('arrowsize', 'The size of the arrowhead'),
             ('hue', 'The color given as a hue.'),
             ('rgbcolor', 'The color as an rgb tuple.'),
             ('width', 'The width of the shaft of the arrow, in points.')]
        """
        return {'width':'The width of the shaft of the arrow, in points.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'arrowshorten':'The length in points to shorten the arrow.',
                'arrowsize':'The size of the arrowhead'}

    def _plot3d_options(self, options=None):
        if options == None:
            options = self.options()
        options = dict(self.options())
        options_3d = {}
        if 'width' in options:
            options_3d['thickness'] = options['width']
            del options['width']
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        EXAMPLE:
            sage: arrow((0,0),(1,1)).plot3d()
        """
        from sage.plot.plot3d.shapes2 import line3d
        options = self._plot3d_options()
        options.update(kwds)
        return line3d([(self.xtail, self.ytail, 0), (self.xhead, self.yhead, 0)], arrow_head=True, **options)

    def _repr_(self):
        """
        Text representation of an arrow graphics primitive.

        EXAMPLES:
            sage: from sage.plot.arrow import Arrow
            sage: Arrow(0,0,2,3,{})._repr_()
            'Arrow from (0.0,0.0) to (2.0,3.0)'
        """
        return "Arrow from (%s,%s) to (%s,%s)"%(self.xtail, self.ytail, self.xhead, self.yhead)

    def _render_on_subplot(self, subplot):
        """
        Render this arrow in a subplot.  This is the key function that
        defines how this arrow graphics primitive is rendered in
        matplotlib's library.

        EXAMPLES:
        This function implicitly ends up rendering this arrow on a matplotlib subplot:
            sage: arrow((0,1), (2,-1))
        """
        options = self.options()
        width = float(options['width'])
        arrowshorten = float(options.get('arrowshorten',0))
        arrowsize = float(options.get('arrowsize',10))
        from matplotlib.arrow_line import ArrowLine
        p = ArrowLine([self.xtail, self.xhead], [self.ytail, self.yhead],  lw=width, arrow='>', arrowsize=arrowsize, arrowshorten=arrowshorten)


        c = to_mpl_color(options['rgbcolor'])
        p._arrowedgecolor=(c)
        p._arrowfacecolor=(c)
        p.set_color(c)
        p.set_solid_capstyle('butt')
        p.set_solid_joinstyle('bevel')
        subplot.add_line(p)


@options(width=2, rgbcolor=(0,0,1))
def arrow(tailpoint, headpoint, **options):
    """
    An arrow from (xmin, ymin) to (xmax, ymax).

    INPUT
        width -- (default 2) the width of the arrow shaft, in points
        rgbcolor -- (default (0,0,1)) the color of the arrow (as an rgb tuple)
        hue -- the color of the arrow (as a number)
        arrowsize -- the size of the arrowhead
        arrowshorten -- the length in points to shorten the arrow

    EXAMPLES:

    A straight, blue arrow
       sage: arrow((1, 1), (3, 3))

    Make a red arrow:
       sage: arrow((-1, -1), (2, 3), rgbcolor=(1,0,0))

    You can change the width of an arrow:
        sage: arrow((1, 1), (3, 3), width=5, arrowsize=15)

    A pretty circle of arrows:
        sage: sum([arrow((0,0), (cos(x),sin(x)), hue=x/(2*pi)) for x in [0..2*pi,step=0.1]]).show(aspect_ratio=1)

    If we want to draw the arrow between objects, for example, the
    boundaries of two lines, we can use the arrowshorten option
    to make the arrow shorter by a certain number of points.
        sage: line([(0,0),(1,0)],thickness=10)+line([(0,1),(1,1)], thickness=10)+arrow((0.5,0),(0.5,1), arrowshorten=10,rgbcolor=(1,0,0))


    """
    from sage.plot.plot import Graphics
    xtail, ytail = tailpoint
    xhead, yhead = headpoint
    g = Graphics()
    g.add_primitive(Arrow(xtail, ytail, xhead, yhead, options=options))
    return g
