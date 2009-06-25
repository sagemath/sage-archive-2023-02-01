"""
Line Plots
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
from sage.plot.primitive import GraphicPrimitive_xydata
from sage.plot.misc import options, rename_keyword
from sage.plot.colors import to_mpl_color

class Line(GraphicPrimitive_xydata):
    """
    Primitive class that initializes the line graphics type.

    EXAMPLES::

        sage: from sage.plot.line import Line
        sage: Line([1,2,7], [1,5,-1], {})
        Line defined by 3 points
    """
    def __init__(self, xdata, ydata, options):
        """
        Initialize a line graphics primitive.

        EXAMPLES::

            sage: from sage.plot.line import Line
            sage: Line([-1,2], [17,4], {'thickness':2})
            Line defined by 2 points
        """
        valid_options = self._allowed_options().keys()
        for opt in options.iterkeys():
            if opt not in valid_options:
                raise RuntimeError("Error in line(): option '%s' not valid."%opt)
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive_xydata.__init__(self, options)

    def _allowed_options(self):
        """
        Displayed the list of allowed line options.

        EXAMPLES::

            sage: from sage.plot.line import Line
            sage: list(sorted(Line([-1,2], [17,4], {})._allowed_options().iteritems()))
            [('alpha', 'How transparent the line is.'),
             ('hue', 'The color given as a hue.'),
             ('linestyle',
              "The style of the line, which is one of '--' (dashed), '-.' (dash dot), '-' (solid), 'steps', ':' (dotted)."),
             ('marker',
              "'0' (tickleft), '1' (tickright), '2' (tickup), '3' (tickdown), '' (nothing), ' ' (nothing), '+' (plus), ',' (pixel), '.' (point), '1' (tri_down), '3' (tri_left), '2' (tri_up), '4' (tri_right), '<' (triangle_left), '>' (triangle_right), 'None' (nothing), 'D' (diamond), 'H' (hexagon2), '_' (hline), '^' (triangle_up), 'd' (thin_diamond), 'h' (hexagon1), 'o' (circle), 'p' (pentagon), 's' (square), 'v' (triangle_down), 'x' (x), '|' (vline)"),
             ('markeredgecolor', 'the markerfacecolor can be any color arg'),
             ('markeredgewidth', 'the size of the markter edge in points'),
             ('markersize', 'the size of the marker in points'),
             ('rgbcolor', 'The color as an rgb tuple.'),
             ('thickness', 'How thick the line is.'),
             ('zorder', 'The layer level in which to draw')]
        """
        return {'alpha':'How transparent the line is.',
                'thickness':'How thick the line is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'linestyle':"The style of the line, which is one of '--' (dashed), '-.' (dash dot), '-' (solid), 'steps', ':' (dotted).",
                'marker':"'0' (tickleft), '1' (tickright), '2' (tickup), '3' (tickdown), '' (nothing), ' ' (nothing), '+' (plus), ',' (pixel), '.' (point), '1' (tri_down), '3' (tri_left), '2' (tri_up), '4' (tri_right), '<' (triangle_left), '>' (triangle_right), 'None' (nothing), 'D' (diamond), 'H' (hexagon2), '_' (hline), '^' (triangle_up), 'd' (thin_diamond), 'h' (hexagon1), 'o' (circle), 'p' (pentagon), 's' (square), 'v' (triangle_down), 'x' (x), '|' (vline)",
                'markersize':'the size of the marker in points',
                'markeredgecolor':'the markerfacecolor can be any color arg',
                'markeredgewidth':'the size of the markter edge in points',
                'zorder':'The layer level in which to draw'
                }

    def _plot3d_options(self, options=None):
        if options == None:
            options = dict(self.options())
        options_3d = {}
        if 'thickness' in options:
            options_3d['thickness'] = options['thickness']
            del options['thickness']
        if 'linestyle' in options:
            if options['linestyle'] != '--':
                raise NotImplementedError, "Invalid 3d line style: %s" % options['linestyle']
            del options['linestyle']
        options_3d.update(GraphicPrimitive_xydata._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        EXAMPLES::

            sage: EllipticCurve('37a').plot(thickness=5).plot3d()
        """
        from sage.plot.plot3d.shapes2 import line3d
        options = self._plot3d_options()
        options.update(kwds)
        return line3d([(x, y, 0) for x, y in zip(self.xdata, self.ydata)], **options)

    def _repr_(self):
        """
        String representation of a line primitive.

        EXAMPLES::

            sage: from sage.plot.line import Line
            sage: Line([-1,2,3,3], [17,4,0,2], {})._repr_()
            'Line defined by 4 points'
        """
        return "Line defined by %s points"%len(self)

    def __getitem__(self, i):
        """
        Extract the i-th element of the line (which is stored as a list of points).

        INPUT:

        - ``i`` -- an integer between 0 and the number of points minus 1

        OUTPUT:

        A 2-tuple of floats.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: line_primitive = L[0]; line_primitive
            Line defined by 4 points
            sage: line_primitive[0]
            (1.0, 2.0)
            sage: line_primitive[2]
            (2.0, 5.0)
            sage: list(line_primitive)
            [(1.0, 2.0), (3.0, -4.0), (2.0, 5.0), (1.0, 2.0)]
        """
        return self.xdata[i], self.ydata[i]

    def __setitem__(self, i, point):
        """
        Set the i-th element of this line (really a sequence of lines
        through given points).

        INPUT:

        - ``i`` -- an integer between 0 and the number of points on the
          line minus 1

        - ``point`` -- a 2-tuple of floats

        EXAMPLES:

        We create a line graphics object $L$ and get ahold of the
        corresponding line graphics primitive::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: line_primitive = L[0]; line_primitive
            Line defined by 4 points

        We then set the 0th point to `(0,0)` instead of `(1,2)`::

            sage: line_primitive[0] = (0,0)
            sage: line_primitive[0]
            (0.0, 0.0)

        Plotting we visibly see the change -- now the line starts at `(0,0)`::

            sage: L
        """
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        r"""
        Return the number of points on this line (where a line is really a sequence
        of line segments through a given list of points).

        EXAMPLES:

        We create a line, then grab the line primitive as \code{L[0]} and compute
        its length::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: len(L[0])
            4
        """
        return len(self.xdata)

    def _render_on_subplot(self, subplot):
        """
        Render this line on a matplotlib subplot.

        INPUT:

        - ``subplot`` -- a matplotlib subplot

        EXAMPLES:

        This implicitly calls this function::

            sage: line([(1,2), (3,-4), (2, 5), (1,2)])
        """
        import matplotlib.lines as lines
        options = dict(self.options())
        del options['alpha']
        del options['thickness']
        del options['rgbcolor']
        p = lines.Line2D(self.xdata, self.ydata, **options)
        options = self.options()
        a = float(options['alpha'])
        p.set_alpha(a)
        p.set_linewidth(float(options['thickness']))
        p.set_color(to_mpl_color(options['rgbcolor']))
        subplot.add_line(p)

def line(points, **kwds):
    """
    Returns either a 2-dimensional or 3-dimensional line depending
    on value of points.

    For information regarding additional arguments, see either line2d?
    or line3d?.

    EXAMPLES::

        sage: line([(0,0), (1,1)])
        sage: line([(0,0,1), (1,1,1)])
    """
    try:
        return line2d(points, **kwds)
    except ValueError:
        from sage.plot.plot3d.shapes2 import line3d
        return line3d(points, **kwds)


@options(alpha=1, rgbcolor=(0,0,1), thickness=1)
def line2d(points, **options):
    r"""
    Create the line through the given list of points.

    Type \code{line2d.options} for a dictionary of the default options for
    lines.  You can change this to change the defaults for all future
    lines.  Use \code{line2d.reset()} to reset to the default options.

    INPUT:

    - ``alpha`` -- How transparent the line is

    - ``thickness`` -- How thick the line is

    - ``rgbcolor`` -- The color as an rgb tuple

    - ``hue`` -- The color given as a hue

    Any MATPLOTLIB line option may also be passed in.  E.g.,

    - ``linestyle`` -- The style of the line, which is one of:

      - ``'--'`` (dashed)
      - ``'-.'`` (dash dot)
      - ``'-'`` (solid)
      - ``'steps'``
      - ``':'`` (dotted)

    - ``marker`` -- The style of the marker, which is one of:

      - ``'0'`` (tickleft)
      - ``'1'`` (tickright)
      - ``'2'`` (tickup)
      - ``'3'`` (tickdown)
      - ``''`` (nothing)
      - ``' '`` (nothing)
      - ``'+'`` (plus)
      - ``','`` (pixel)
      - ``'.'`` (point)
      - ``'1'`` (tri_down)
      - ``'3'`` (tri_left)
      - ``'2'`` (tri_up)
      - ``'4'`` (tri_right)
      - ``'<'`` (triangle_left)
      - ``'>'`` (triangle_right)
      - ``'None'`` (nothing)
      - ``'D'`` (diamond)
      - ``'H'`` (hexagon2)
      - ``'_'`` (hline)
      - ``'\^'`` (triangle_up)
      - ``'d'`` (thin_diamond)
      - ``'h'`` (hexagon1)
      - ``'o'`` (circle)
      - ``'p'`` (pentagon)
      - ``'s'`` (square)
      - ``'v'`` (triangle_down)
      - ``'x'`` (x)
      - ``'|'`` (vline)

    - ``markersize`` -- the size of the marker in points

    - ``markeredgecolor`` -- the markerfacecolor can be any color arg

    - ``markeredgewidth`` -- the size of the markter edge in points

    EXAMPLES:

    A blue conchoid of Nicomedes::

        sage: L = [[1+5*cos(pi/2+pi*i/100), tan(pi/2+pi*i/100)*(1+5*cos(pi/2+pi*i/100))] for i in range(1,100)]
        sage: line(L, rgbcolor=(1/4,1/8,3/4))

    A line with 2 complex points::

        sage: i = CC.0
        sage: line([1+i, 2+3*i])

    A blue hypotrochoid (3 leaves)::

        sage: n = 4; h = 3; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: line(L, rgbcolor=(1/4,1/4,3/4))

    A blue hypotrochoid (4 leaves)::

        sage: n = 6; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: line(L, rgbcolor=(1/4,1/4,3/4))

    A red limacon of Pascal::

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,101)]
        sage: line(L, rgbcolor=(1,1/4,1/2))

    A light green trisectrix of Maclaurin::

        sage: L = [[2*(1-4*cos(-pi/2+pi*i/100)^2),10*tan(-pi/2+pi*i/100)*(1-4*cos(-pi/2+pi*i/100)^2)] for i in range(1,100)]
        sage: line(L, rgbcolor=(1/4,1,1/8))

    A green lemniscate of Bernoulli::

        sage: cosines = [cos(-pi/2+pi*i/100) for i in range(201)]
        sage: v = [(1/c, tan(-pi/2+pi*i/100)) for i,c in enumerate(cosines) if c != 0]
        sage: L = [(a/(a^2+b^2), b/(a^2+b^2)) for a,b in v]
        sage: line(L, rgbcolor=(1/4,3/4,1/8))

    A red plot of the Jacobi elliptic function `\text{sn}(x,2)`, `-3 < x < 3`::

        sage: L = [(i/100.0, jacobi('sn', i/100.0 ,2.0)) for i in range(-300,300,30)]
        sage: line(L, rgbcolor=(3/4,1/4,1/8))

    A red plot of `J`-Bessel function `J_2(x)`, `0 < x < 10`::

        sage: L = [(i/10.0, bessel_J(2,i/10.0)) for i in range(100)]
        sage: line(L, rgbcolor=(3/4,1/4,5/8))


    A purple plot of the Riemann zeta function `\zeta(1/2 + it)`, `0 < t < 30`::

        sage: i = CDF.gen()
        sage: v = [zeta(0.5 + n/10 * i) for n in range(300)]
        sage: L = [(z.real(), z.imag()) for z in v]
        sage: line(L, rgbcolor=(3/4,1/2,5/8))

    A purple plot of the Hasse-Weil `L`-function `L(E, 1 + it)`, `-1 < t < 10`::

        sage: E = EllipticCurve('37a')
        sage: vals = E.lseries().values_along_line(1-I, 1+10*I, 100) # critical line
        sage: L = [(z[1].real(), z[1].imag()) for z in vals]
        sage: line(L, rgbcolor=(3/4,1/2,5/8))

    A red, blue, and green "cool cat"::

        sage: G = plot(-cos(x), -2, 2, thickness=5, rgbcolor=(0.5,1,0.5))
        sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))
        sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
        sage: G + P + Q   # show the plot

    A line with no points or one point::

        sage: line([])
        sage: line([(1,1)])

    """
    from sage.plot.plot import Graphics, xydata_from_point_list
    xdata, ydata = xydata_from_point_list(points)
    g = Graphics()
    g.add_primitive(Line(xdata, ydata, options))
    return g
