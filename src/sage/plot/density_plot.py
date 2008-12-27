#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#                          Arnaud Bergeron <abergeron@gmail.com>
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
from sage.plot.misc import options, rename_keyword, to_mpl_color, get_cmap
from sage.misc.misc import verbose, xsrange

class DensityPlot(GraphicPrimitive):
    """
    Primitive class that initializes the
    density_plot graphics type
    """
    def __init__(self, xy_data_array, xrange, yrange, options):
        """
        TESTS:
            sage: dp = density_plot(x, (-2,3), (1,10))
        """
        self.xrange = xrange
        self.yrange = yrange
        self.xy_data_array = xy_data_array
        self.xy_array_row = len(xy_data_array)
        self.xy_array_col = len(xy_data_array[0])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES:
            sage: x,y = var('x,y')
            sage: f = x^2 + y^2
            sage: d = density_plot(f, (3, 6), (3, 6))[0].get_minmax_data()
            sage: d['xmin']
            3.0
            sage: d['ymin']
            3.0

        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xrange, self.yrange, dict=True)

    def _allowed_options(self):
        """
        TESTS:
            sage: isinstance(density_plot(x, (-2,3), (1,10))[0]._allowed_options(), dict)
            True
        """
        return {'plot_points':'How many points to use for plotting precision',
                'cmap':"A colormap (type cmap_help() for more information).",
                'interpolation':'What interpolation method to use'}

    def _repr_(self):
        """
        TESTS:
            sage: isinstance(density_plot(x, (-2,3), (1,10))[0]._repr_(), str)
            True
        """
        return "DensityPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        options = self.options()
        cmap = get_cmap(options['cmap'])

        x0,x1 = float(self.xrange[0]), float(self.xrange[1])
        y0,y1 = float(self.yrange[0]), float(self.yrange[1])

        subplot.imshow(self.xy_data_array, origin='lower', cmap=cmap, extent=(x0,x1,y0,y1), interpolation=options['interpolation'])

@options(plot_points=25, cmap='gray', interpolation='catrom')
def density_plot(f, xrange, yrange, **options):
    r"""

    \code{density_plot} takes a function of two variables, $f(x,y)$
    and plots the height of of the function over the specified
    xrange and yrange as demonstrated below.

      density_plot(f, (xmin, xmax), (ymin, ymax), ...)

    INPUT:
        f -- a function of two variables
        (xmin, xmax) -- 2-tuple, the range of x values OR 3-tuple (x,xmin,xmax)
        (ymin, ymax) -- 2-tuple, the range of y values OR 3-tuple (y,ymin,ymax)
    The following inputs must all be passed in as named parameters:
        plot_points   -- integer (default: 25); number of points to plot
                         in each direction of the grid
        cmap          -- a colormap (type cmap_help() for more information).
        interpolation -- string (default: 'catrom'), the interpolation
                         method to use: bilinear, bicubic, spline16, spline36,
                         quadric, gaussian, sinc, bessel, mitchell, lanczos,
                         catrom, hermite, hanning, hamming, kaiser


    EXAMPLES:
    Here we plot a simple function of two variables:
        sage: x,y = var('x,y')
        sage: density_plot(sin(x)*sin(y), (-2, 2), (-2, 2))


    Here we change the ranges and add some options:
        sage: density_plot((x^2)*cos(x*y), (-10, 5), (-5, 5), interpolation='sinc', plot_points=100)


    An even more complicated plot.
        sage: density_plot(sin(x^2 + y^2)*cos(x)*sin(y), (-4, 4), (-4, 4), cmap='jet', plot_points=100)


    Some elliptic curves, but with symbolic endpoints.  In the first
    example, the plot is rotated 90 degrees because we switch the
    variables x,y.
        sage: density_plot(y^2 + 1 - x^3 - x, (y,-pi,pi), (x,-pi,pi))
        sage: density_plot(y^2 + 1 - x^3 - x, (-pi,pi), (-pi,pi))
    """
    from sage.plot.plot import Graphics, setup_for_eval_on_grid
    g, xstep, ystep, xrange, yrange = setup_for_eval_on_grid([f], xrange, yrange, options['plot_points'])
    g = g[0]
    xy_data_array = [[g(x, y) for x in xsrange(xrange[0], xrange[1], xstep)]
                              for y in xsrange(yrange[0], yrange[1], ystep)]

    g = Graphics()
    g.add_primitive(DensityPlot(xy_data_array, xrange, yrange, options))
    return g
