"""
Streamline Plots
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
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options
from sage.arith.srange import xsrange

class StreamlinePlot(GraphicPrimitive):
    """
    Primitive class that initializes the StreamlinePlot graphics type
    """
    def __init__(self, xpos_array, ypos_array, xvec_array, yvec_array, options):
        """
        Create the graphics primitive StreamlinePlot. This sets options
        and the array to be plotted as attributes.

        EXAMPLES::

            sage: x, y = var('x y')
            sage: R = streamline_plot((sin(x), cos(y)), (x,0,1), (y,0,1), plot_points=2)
            sage: r = R[0]
            sage: r.options()['plot_points']
            2
            sage: r.xpos_array
            array([ 0.,  1.])
            sage: r.yvec_array
            masked_array(data =
             [[1.0 1.0]
             [0.5403023058681398 0.5403023058681398]],
                         mask =
             [[False False]
             [False False]],
                   fill_value = 1e+20)
            <BLANKLINE>

        TESTS:

        We test dumping and loading a plot::

            sage: x, y = var('x y')
            sage: P = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3))
            sage: Q = loads(dumps(P))

        """
        self.xpos_array = xpos_array
        self.ypos_array = ypos_array
        self.xvec_array = xvec_array
        self.yvec_array = yvec_array
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: x, y = var('x y')
            sage: d = streamline_plot((.01*x, x+y), (x,10,20), (y,10,20))[0].get_minmax_data()
            sage: d['xmin']
            10.0
            sage: d['ymin']
            10.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xpos_array, self.ypos_array, dict=True)

    def _allowed_options(self):
        """
        Returns a dictionary with allowed options for StreamlinePlot.

        EXAMPLES::

            sage: x, y = var('x y')
            sage: P = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3))
            sage: d = P[0]._allowed_options()
            sage: d['pivot']
            'Where the arrow should be placed in relation to the point (tail, middle, tip)'
        """
        return {'plot_points': 'How many points to use for plotting precision',
                'color': 'The color of the arrows',
                'zorder': 'The layer level in which to draw'}

    def _repr_(self):
        """
        String representation of StreamlinePlot graphics primitive.

        EXAMPLES::

            sage: x, y = var('x y')
            sage: P = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3))
            sage: P[0]
            StreamlinePlot defined by a 20 x 20 vector grid

        TESTS:

            sage: x, y = var('x y')
            sage: P = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3), wrong_option='nonsense')
            sage: P[0].options()['plot_points']
            verbose 0 (...: primitive.py, options) WARNING: Ignoring option 'wrong_option'=nonsense
            verbose 0 (...: primitive.py, options)
            The allowed options for StreamlinePlot defined by a 20 x 20 vector grid are:
                color          The color of the arrows
                plot_points    How many points to use for plotting precision
                zorder         The layer level in which to draw
            <BLANKLINE>
            20

        """
        return "StreamlinePlot defined by a {} x {} vector grid".format(
               self._options['plot_points'], self._options['plot_points'])

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: x, y = var('x y')
            sage: P = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3))
        """
        options = self.options()
        streamplot_options = options.copy()
        streamplot_options.pop('plot_points')
        subplot.streamplot(self.xpos_array, self.ypos_array,
                       self.xvec_array, self.yvec_array,
                       **streamplot_options)


@options(plot_points=20, frame=True)
def streamline_plot(f_g, xrange, yrange, **options):
    r"""
    ``streamline_plot`` takes two functions of two variables ``xvar`` and ``yvar``
    (for instance, if the variables are `x` and `y`, take `(f(x,y), g(x,y))`)
    and plots vector steamlines of the function over the specified ranges, with
    ``xrange`` being of ``xvar`` between ``xmin`` and ``xmax``, and ``yrange`` similarly
    (see below).

    ``streamline_plot((f, g), (xvar,xmin,xmax), (yvar,ymin,ymax))``

    EXAMPLES:

    Plot some vector fields involving sin and cos::

        sage: x, y = var('x y')
        sage: streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var('x y')
        g = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3))
        sphinx_plot(g)

    ::

        sage: streamline_plot((y, (cos(x)-2) * sin(x)), (x,-pi,pi), (y,-pi,pi))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var('x y')
        g = streamline_plot((y, (cos(x)-2) * sin(x)), (x,-pi,pi), (y,-pi,pi))
        sphinx_plot(g)

    We ignore function values that are infinite or NaN::

        sage: x, y = var('x y')
        sage: streamline_plot((-x/sqrt(x^2+y^2), -y/sqrt(x^2+y^2)), (x,-10,10), (y,-10,10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var('x y')
        g = streamline_plot((-x/sqrt(x**2+y**2), -y/sqrt(x**2+y**2)), (x,-10,10), (y,-10,10))
        sphinx_plot(g)

    Extra options will get passed on to show(), as long as they are valid::

        sage: streamline_plot((x, y), (x,-2,2), (y,-2,2), xmax=10)
        Graphics object consisting of 1 graphics primitive
        sage: streamline_plot((x, y), (x,-2,2), (y,-2,2)).show(xmax=10) # These are equivalent

    .. PLOT::

        x, y = var('x y')
        g = streamline_plot((x, y), (x,-2,2), (y,-2,2), xmax=10)
        sphinx_plot(g)

    """
    (f,g) = f_g
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    z, ranges = setup_for_eval_on_grid([f,g], [xrange,yrange], options['plot_points'])
    f, g = z

    xpos_array, ypos_array, xvec_array, yvec_array = [], [], [], []
    for x in xsrange(*ranges[0], include_endpoint=True):
        xpos_array.append(x)
    for y in xsrange(*ranges[1], include_endpoint=True):
        ypos_array.append(y)
        xvec_row, yvec_row = [], []
        for x in xsrange(*ranges[0], include_endpoint=True):
            xvec_row.append(f(x, y))
            yvec_row.append(g(x, y))
        xvec_array.append(xvec_row)
        yvec_array.append(yvec_row)

    import numpy
    xpos_array = numpy.array(xpos_array, dtype=float)
    ypos_array = numpy.array(ypos_array, dtype=float)
    xvec_array = numpy.ma.masked_invalid(numpy.array(xvec_array, dtype=float))
    yvec_array = numpy.ma.masked_invalid(numpy.array(yvec_array, dtype=float))
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(StreamlinePlot(xpos_array, ypos_array,
                                   xvec_array, yvec_array, options))
    return g
