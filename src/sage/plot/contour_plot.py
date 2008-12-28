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
from sage.plot.misc import options, rename_keyword, to_mpl_color, get_cmap
from sage.misc.misc import verbose, xsrange

class ContourPlot(GraphicPrimitive):
    """
    Primitive class that initializes the
    contour_plot graphics type
    """
    def __init__(self, xy_data_array, xrange, yrange, options):
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
            sage: d = contour_plot(f, (3, 6), (3, 6))[0].get_minmax_data()
            sage: d['xmin']
            3.0
            sage: d['ymin']
            3.0

        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xrange, self.yrange, dict=True)

    def _allowed_options(self):
        return {'plot_points':'How many points to use for plotting precision',
                'cmap':"""the name of a predefined colormap,
                        a list of colors or an instance of a
                        matplotlib Colormap.""",
                       'fill':'Fill contours or not',
                'contours':"""Either an integer specifying the number of
                       contour levels, or a sequence of numbers giving
                       the actual contours to use."""}

    def _repr_(self):
        return "ContourPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        from sage.rings.integer import Integer
        options = self.options()
        fill = options['fill']
        contours = options['contours']
        cmap = get_cmap(options['cmap'])

        x0,x1 = float(self.xrange[0]), float(self.xrange[1])
        y0,y1 = float(self.yrange[0]), float(self.yrange[1])
        if fill:
            if contours is None:
                subplot.contourf(self.xy_data_array, cmap=cmap, extent=(x0,x1,y0,y1))
            elif isinstance(contours, (int, Integer)):
                subplot.contourf(self.xy_data_array, int(contours), cmap=cmap, extent=(x0,x1,y0,y1))
            else:
                subplot.contourf(self.xy_data_array, contours, cmap=cmap, extent=(x0,x1,y0,y1))
        else:
            if contours is None:
                subplot.contour(self.xy_data_array, cmap=cmap, extent=(x0,x1,y0,y1))
            elif isinstance(contours, (int, Integer)):
                subplot.contour(self.xy_data_array, int(contours), cmap=cmap, extent=(x0,x1,y0,y1))
            else:
                subplot.contour(self.xy_data_array, contours, cmap=cmap, extent=(x0,x1,y0,y1))

@options(plot_points=25, fill=True, cmap='gray', contours=None)
def contour_plot(f, xrange, yrange, **options):
    r"""

    \code{contour_plot} takes a function of two variables, $f(x,y)$
    and plots contour lines of the function over the specified
    xrange and yrange as demonstrated below.

      contour_plot(f, (xmin, xmax), (ymin, ymax), ...)

    INPUT:
        f -- a function of two variables
        (xmin, xmax) -- 2-tuple, the range of x values OR 3-tuple (x,xmin,xmax)
        (ymin, ymax) -- 2-tuple, the range of y values OR 3-tuple (y,ymin,ymax)
    The following inputs must all be passed in as named parameters:
        plot_points  -- integer (default: 25); number of points to plot
                        in each direction of the grid
        fill         -- bool (default: True), whether to color in the area
                        between contour lines
        cmap         -- a colormap (default: 'gray'), the name of
                        a predefined colormap, a list of colors
                        or an instance of a matplotlib Colormap.
                        Type: import matplotlib.cm; matplotlib.cm.datad.keys()
                        for available colormap names.
        contours     -- integer or list of numbers (default: None):
                        If a list of numbers is given, then this specifies
                        the contour levels to use.  If an integer is given,
                        then this many contour lines are used, but the
                        exact levels are determined automatically.
                        If None is passed (or the option is not given),
                        then the number of contour lines is determined
                        automatically, and is usually about 5.


    EXAMPLES:

    Here we plot a simple function of two variables:
        sage: x,y = var('x,y')
        sage: contour_plot(cos(x^2+y^2), (-4, 4), (-4, 4))


    Here we change the ranges and add some options:
        sage: contour_plot((x^2)*cos(x*y), (-10, 5), (-5, 5), fill=False, plot_points=100)


    An even more complicated plot.
        sage: contour_plot(sin(x^2 + y^2)*cos(x)*sin(y), (-4, 4), (-4, 4),plot_points=100)

    Some elliptic curves, but with symbolic endpoints.  In the first
    example, the plot is rotated 90 degrees because we switch the
    variables x,y.
        sage: contour_plot(y^2 + 1 - x^3 - x, (y,-pi,pi), (x,-pi,pi))
        sage: contour_plot(y^2 + 1 - x^3 - x, (-pi,pi), (-pi,pi))


    We can play with the contour levels.
        sage: f = x^2 + y^2
        sage: contour_plot(f, (-2, 2), (-2, 2))
        sage: contour_plot(f, (-2, 2), (-2, 2), contours=2, cmap=[(1,0,0), (0,1,0), (0,0,1)])
        sage: contour_plot(f, (-2, 2), (-2, 2), contours=(0.1, 1.0, 1.2, 1.4), cmap='hsv')
        sage: contour_plot(f, (-2, 2), (-2, 2), contours=(1.0,), fill=False)

    """
    from sage.plot.plot import Graphics, setup_for_eval_on_grid
    g, xstep, ystep, xrange, yrange = setup_for_eval_on_grid([f], xrange, yrange, options['plot_points'])
    g = g[0]
    xy_data_array = [[g(x, y) for x in xsrange(xrange[0], xrange[1], xstep)]
                              for y in xsrange(yrange[0], yrange[1], ystep)]

    g = Graphics()
    g.add_primitive(ContourPlot(xy_data_array, xrange, yrange, options))
    return g

@options(contours=(0,0), fill=False)
def implicit_plot(f, xrange, yrange, **options):
    r"""
    \code{implicit_plot} takes a function of two variables, $f(x,y)$
    and plots the curve $f(x,y)=0$ over the specified
    xrange and yrange as demonstrated below.

      implicit_plot(f, (xmin, xmax), (ymin, ymax), ...)

    INPUT:
        f -- a function of two variables
        (xmin, xmax) -- 2-tuple, the range of x values
        (ymin, ymax) -- 2-tuple, the range of y values
    The following inputs must all be passed in as named parameters:
        plot_points  -- integer (default: 25); number of points to plot
                        in each direction of the grid
        fill         -- boolean (default: False); if True, fill the region $f(x,y)<0$.


    EXAMPLES:

    A simple circle with a radius of 2:
        sage: var("x y")
        (x, y)
        sage: implicit_plot(x^2+y^2-2, (-3,3), (-3,3)).show(aspect_ratio=1)

    We can define a level-$n$ approximation of the boundary of the
    Mandelbrot set.
        sage: def mandel(n):
        ...       c = polygen(CDF, 'c')
        ...       z = 0
        ...       for i in range(n):
        ...           z = z*z + c
        ...       def f(x, y):
        ...           val = z(CDF(x, y))
        ...           return val.norm() - 4
        ...       return f

    The first-level approximation is just a circle.
        sage: implicit_plot(mandel(1), (-3, 3), (-3, 3)).show(aspect_ratio=1)

    A third-level approximation starts to get interesting.
        sage: implicit_plot(mandel(3), (-2, 1), (-1.5, 1.5)).show(aspect_ratio=1)

    The seventh-level approximation is a degree 64 polynomial, and
    implicit_plot does a pretty good job on this part of the curve.
    (plot_points=200 looks even better, but it's about 16 times slower.)
        sage: implicit_plot(mandel(7), (-0.3, 0.05), (-1.15, -0.9),plot_points=50).show(aspect_ratio=1)
    """
    return contour_plot(f, xrange, yrange, **options)

@options(plot_points=25, incol='blue', outcol='white', bordercol=None)
def region_plot(f, xrange, yrange, plot_points, incol, outcol, bordercol):
    r"""
    \code{region_plot} takes a boolean function of two variables, $f(x,y)$
    and plots the region where f is True over the specified
    xrange and yrange as demonstrated below.

      region_plot(f, (xmin, xmax), (ymin, ymax), ...)

    INPUT:
        f -- a boolean function of two variables
        (xmin, xmax) -- 2-tuple, the range of x values OR 3-tuple (x,xmin,xmax)
        (ymin, ymax) -- 2-tuple, the range of y values OR 3-tuple (y,ymin,ymax)
        plot_points  -- integer (default: 25); number of points to plot
                        in each direction of the grid
        incol        -- a color (default: 'blue'), the color inside the region
        outcol       -- a color (default: 'white'), the color of the outside of the region
        bordercol    -- a color (default: None), the color of the border (incol if not specified)

    EXAMPLES:

    Here we plot a simple function of two variables:
        sage: x,y = var('x,y')
        sage: region_plot(cos(x^2+y^2)<0, (-3, 3), (-3, 3))

    Here we play with the colors:
        sage: region_plot(x^2+y^3<2, (-2, 2), (-2, 2), incol='lightblue', bordercol='gray')

    An even more complicated plot.
        sage: region_plot(sin(x)*sin(y) >=1/4, (-10,10), (-10,10), incol='yellow', bordercol='black', plot_points=100).show(aspect_ratio=1)
    """
    from matplotlib.colors import ListedColormap
    import operator
    incol = rgbcolor(incol)
    outcol = rgbcolor(outcol)
    cmap = ListedColormap([incol, outcol])
    cmap.set_over(outcol)
    cmap.set_under(incol)
    op = f.operator()
    if op is operator.gt or op is operator.ge:
        g = f.rhs() - f.lhs()
    else:
        g = f.lhs() - f.rhs()

    from sage.plot.plot import Graphics
    res = Graphics()
    if op is not operator.eq:
        res += contour_plot(g, xrange, yrange, plot_points=plot_points,
                           contours=[-1e307, 0, 1e307], cmap=cmap, fill=True)

    if op is operator.ge or op is operator.le or op is operator.eq or bordercol is not None:
        if bordercol is None: bordercol = incol
        res += contour_plot(g, xrange, yrange, plot_points=plot_points, contours=[0], cmap=[bordercol], fill=False)

    return res
