r"""
2D Plotting

Sage provides extensive 2D plotting functionality. The underlying
rendering is done using the matplotlib Python library.

The following graphics primitives are supported:


-  :func:`~sage.plot.arrow.arrow` - an arrow from a min point to a max point.

-  :func:`~sage.plot.circle.circle` - a circle with given radius

-  :func:`~sage.plot.ellipse.ellipse` - an ellipse with given radii
   and angle

-  :func:`~sage.plot.arc.arc` - an arc of a circle or an ellipse

-  :func:`~sage.plot.disk.disk` - a filled disk (i.e. a sector or wedge of a circle)

-  :func:`~sage.plot.line.line` - a line determined by a sequence of points (this need not
   be straight!)

-  :func:`~sage.plot.point.point` - a point

-  :func:`~sage.plot.text.text` - some text

-  :func:`~sage.plot.polygon.polygon` - a filled polygon


The following plotting functions are supported:


-  :func:`plot` - plot of a function or other Sage object (e.g., elliptic
   curve).

-  :func:`parametric_plot`

-  :func:`~sage.plot.contour_plot.implicit_plot`

-  :func:`polar_plot`

-  :func:`~sage.plot.contour_plot.region_plot`

-  :func:`list_plot`

-  :func:`~sage.plot.scatter_plot.scatter_plot`

-  :func:`~sage.plot.bar_chart.bar_chart`

-  :func:`~sage.plot.contour_plot.contour_plot`

-  :func:`~sage.plot.density_plot.density_plot`

-  :func:`~sage.plot.plot_field.plot_vector_field`

-  :func:`~sage.plot.plot_field.plot_slope_field`

-  :func:`~sage.plot.matrix_plot.matrix_plot`

-  :func:`~sage.plot.complex_plot.complex_plot`

-  :func:`graphics_array`


The following miscellaneous Graphics functions are included:


-  :func:`Graphics`

-  :func:`is_Graphics`

-  :func:`~sage.plot.colors.hue`


Type ``?`` after each primitive in Sage for help and examples.

EXAMPLES:

We draw a curve::

    sage: plot(x^2, (x,0,5))

We draw a circle and a curve::

    sage: circle((1,1), 1) + plot(x^2, (x,0,5))

Notice that the aspect ratio of the above plot makes the plot very tall because
the plot adopts the default aspect ratio of the circle (to make the circle appear
like a circle).  We can change the aspect ratio to be what we normally expect for a plot
by explicitly asking for an 'automatic' aspect ratio::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), aspect_ratio='automatic')

The aspect ratio describes the apparently height/width ratio of a unit square.  If you want the vertical units to be twice as big as the horizontal units, specify an aspect ratio of 2::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), aspect_ratio=2)

The ``figsize`` option adjusts the figure size.  The default figsize is 4.  To make a figure that is roughly twice as big, use ``figsize=8``::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), figsize=8)

You can also give separate horizontal and vertical dimensions::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), figsize=[4,8])

Note that the axes will not cross if the data is not on both sides of
both axes, even if it is quite close::

    sage: plot(x^3,(x,1,10))

When the labels have quite different orders of magnitude or are very
large, scientific notation (the `e` notation for powers of ten) is used::

    sage: plot(x^2,(x,480,500))  # no scientific notation

::

    sage: plot(x^2,(x,300,500))  # scientific notation on y-axis

But you can fix your own tick labels, if you know what to expect and
have a preference::

    sage: plot(x^2,(x,300,500),ticks=[None,50000])

We construct a plot involving several graphics objects::

    sage: G = plot(cos(x), (x, -5, 5), thickness=5, color='green')
    sage: P = polygon([[1,2], [5,6], [5,0]], color='red')
    sage: G + P

Next we construct the reflection of the above polygon about the
`y`-axis by iterating over the list of first-coordinates of
the first graphic element of `P` (which is the actual
Polygon; note that `P` is a Graphics object, which consists
of a single polygon)::

    sage: Q = polygon([(-x,y) for x,y in P[0]], color='blue')
    sage: Q   # show it

We combine together different graphics objects using "+"::

    sage: H = G + P + Q
    sage: print H
    Graphics object consisting of 3 graphics primitives
    sage: type(H)
    <class 'sage.plot.plot.Graphics'>
    sage: H[1]
    Polygon defined by 3 points
    sage: list(H[1])
    [(1.0, 2.0), (5.0, 6.0), (5.0, 0.0)]
    sage: H       # show it

We can put text in a graph::

    sage: L = [[cos(pi*i/100)^3,sin(pi*i/100)] for i in range(200)]
    sage: p = line(L, rgbcolor=(1/4,1/8,3/4))
    sage: t = text('A Bulb', (1.5, 0.25))
    sage: x = text('x axis', (1.5,-0.2))
    sage: y = text('y axis', (0.4,0.9))
    sage: g = p+t+x+y
    sage: g.show(xmin=-1.5, xmax=2, ymin=-1, ymax=1)

We plot the Riemann zeta function along the critical line and see
the first few zeros::

    sage: i = CDF.0      # define i this way for maximum speed.
    sage: p1 = plot(lambda t: arg(zeta(0.5+t*i)), 1,27,rgbcolor=(0.8,0,0))
    sage: p2 = plot(lambda t: abs(zeta(0.5+t*i)), 1,27,color=hue(0.7))
    sage: print p1 + p2
    Graphics object consisting of 2 graphics primitives
    sage: p1 + p2    # display it

Many concentric circles shrinking toward the origin::

    sage: show(sum(circle((i,0), i, hue=sin(i/10)) for i in [10,9.9,..,0]))

Here is a pretty graph::

    sage: g = Graphics()
    sage: for i in range(60):
    ...    p = polygon([(i*cos(i),i*sin(i)), (0,i), (i,0)],\
    ...                color=hue(i/40+0.4), alpha=0.2)
    ...    g = g + p
    ...
    sage: g.show(dpi=200, axes=False)

Another graph::

    sage: x = var('x')
    sage: P = plot(sin(x)/x, -4,4, color='blue') + \
    ...       plot(x*cos(x), -4,4, color='red') + \
    ...       plot(tan(x),-4,4, color='green')
    ...
    sage: P.show(ymin=-pi,ymax=pi)

PYX EXAMPLES: These are some examples of plots similar to some of
the plots in the PyX (http://pyx.sourceforge.net) documentation:

Symbolline::

    sage: y(x) = x*sin(x^2)
    sage: v = [(x, y(x)) for x in [-3,-2.95,..,3]]
    sage: show(points(v, rgbcolor=(0.2,0.6, 0.1), pointsize=30) + plot(spline(v), -3.1, 3))

Cycliclink::

    sage: x = var('x')
    sage: g1 = plot(cos(20*x)*exp(-2*x), 0, 1)
    sage: g2 = plot(2*exp(-30*x) - exp(-3*x), 0, 1)
    sage: show(graphics_array([g1, g2], 2, 1), xmin=0)

Pi Axis::

    sage: g1 = plot(sin(x), 0, 2*pi)
    sage: g2 = plot(cos(x), 0, 2*pi, linestyle = "--")
    sage: (g1+g2).show(ticks=pi/6, tick_formatter=pi)  # show their sum, nicely formatted

An illustration of integration::

    sage: f(x) = (x-3)*(x-5)*(x-7)+40
    sage: P = line([(2,0),(2,f(2))], color='black')
    sage: P += line([(8,0),(8,f(8))], color='black')
    sage: P += polygon([(2,0),(2,f(2))] + [(x, f(x)) for x in [2,2.1,..,8]] + [(8,0),(2,0)],  rgbcolor=(0.8,0.8,0.8),aspect_ratio='automatic')
    sage: P += text("$\\int_{a}^b f(x) dx$", (5, 20), fontsize=16, color='black')
    sage: P += plot(f, (1, 8.5), thickness=3)
    sage: P    # show the result



NUMERICAL PLOTTING:

Sage includes Matplotlib, which provides 2D plotting with an interface
that is a likely very familiar to people doing numerical
computation. For example,

::

    sage: from pylab import *
    sage: t = arange(0.0, 2.0, 0.01)
    sage: s = sin(2*pi*t)
    sage: P = plot(t, s, linewidth=1.0)
    sage: xl = xlabel('time (s)')
    sage: yl = ylabel('voltage (mV)')
    sage: t = title('About as simple as it gets, folks')
    sage: grid(True)
    sage: savefig(os.path.join(SAGE_TMP, 'sage.png'))

We test that ``imshow`` works as well, verifying that
Trac ticket 2900 is fixed in Matplotlib.

::

    sage: imshow([[(0,0,0)]])
    <matplotlib.image.AxesImage object at ...>
    sage: savefig(os.path.join(SAGE_TMP, 'foo.png'))

Since the above overwrites many Sage plotting functions, we reset
the state of Sage, so that the examples below work!

::

    sage: reset()

See http://matplotlib.sourceforge.net for complete documentation
about how to use Matplotlib.

TESTS: We test dumping and loading a plot.

::

    sage: p = plot(sin(x), (x, 0,2*pi))
    sage: Q = loads(dumps(p))

Verify that a clean sage startup does *not* import matplotlib::

    sage: os.system("sage -c \"if 'matplotlib' in sys.modules: sys.exit(1)\"")
    0

AUTHORS:

- Alex Clemesha and William Stein (2006-04-10): initial version

- David Joyner: examples

- Alex Clemesha (2006-05-04) major update

- William Stein (2006-05-29): fine tuning, bug fixes, better server
  integration

- William Stein (2006-07-01): misc polish

- Alex Clemesha (2006-09-29): added contour_plot, frame axes, misc
  polishing

- Robert Miller (2006-10-30): tuning, NetworkX primitive

- Alex Clemesha (2006-11-25): added plot_vector_field, matrix_plot,
  arrow, bar_chart, Axes class usage (see axes.py)

- Bobby Moretti and William Stein (2008-01): Change plot to specify
  ranges using the (varname, min, max) notation.

- William Stein (2008-01-19): raised the documentation coverage from a
  miserable 12 percent to a 'wopping' 35 percent, and fixed and
  clarified numerous small issues.

- Jason Grout (2009-09-05): shifted axes and grid functionality over
  to matplotlib; fixed a number of smaller issues.

- Jason Grout (2010-10): rewrote aspect ratio portions of the code

"""

############################################################################
#  Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com> and William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

import os, types

from sage.structure.sage_object import SageObject

## IMPORTANT: Do *not* import matplotlib at module scope.  It takes a
## surprisingly long time to initialize itself.  It's better if it is
## imported in functions, so it only gets started if it is actually
## going to be used.

ALLOWED_EXTENSIONS = ['.eps', '.pdf', '.png', '.ps', '.sobj', '.svg']
#DEFAULT_FIGSIZE=(6, 3.70820393249937)
DEFAULT_DPI = 100
EMBEDDED_MODE = False
DOCTEST_MODE = False
import sage.misc.misc
from sage.misc.misc import srange
DOCTEST_MODE_FILE = os.path.join(sage.misc.misc.SAGE_TMP, 'test.png')
SHOW_DEFAULT = True

def show_default(default=None):
    r"""
    Set the default for showing plots using any plot commands. If
    called with no arguments, returns the current default.

    If this is ``True`` (the default) then any plot object
    when displayed will be displayed as an actual plot instead of text,
    i.e., the show command is not needed.

    EXAMPLES: The default starts out as ``True``::

        sage: show_default()
        True

    We set it to ``False``.

    ::

        sage: show_default(False)

    We see that it is ``False``.

    ::

        sage: show_default()
        False

    Now plot commands will not display their plots by default.

    Turn back on default display.

    ::

        sage: show_default(True)
    """
    global SHOW_DEFAULT
    if default is None:
        return SHOW_DEFAULT
    SHOW_DEFAULT = bool(default)

# If do_verify is True, options are checked when drawing a
# GraphicsPrimitive.  See primitive.py
do_verify = True

from sage.misc.randstate import current_randstate #for plot adaptive refinement
import os #for viewing and writing images
from math import sin, cos, pi #for polar_plot
from sage.structure.sage_object import SageObject

from sage.ext.fast_eval import fast_float, fast_float_constant, is_fast_float

from sage.misc.html import html

from sage.misc.decorators import options, suboptions, rename_keyword

from colors import hue, rainbow, rgbcolor, Color, to_mpl_color

import operator

############### WARNING ###
# Try not to import any matplotlib stuff here -- matplotlib is
# slow to import.  (I did benchmarking and found that by not
# importing here, and instead importing when needed below, that
# Sage startup times are much improved.)  - William
###############

def is_Graphics(x):
    """
    Return True if `x` is a Graphics object.

    EXAMPLES::

        sage: from sage.plot.plot import is_Graphics
        sage: is_Graphics(1)
        False
        sage: is_Graphics(disk((0.0, 0.0), 1, (0, pi/2)))
        True
    """
    return isinstance(x, Graphics)

class Graphics(SageObject):
    """
    The Graphics object is an empty list of graphics objects It is
    useful to use this object when initializing a for loop where
    different graphics object will be added to the empty object.

    EXAMPLES::

        sage: G = Graphics(); print G
        Graphics object consisting of 0 graphics primitives
        sage: c = circle((1,1), 1)
        sage: G+=c; print G
        Graphics object consisting of 1 graphics primitive

    Here we make a graphic of embedded isosceles triangles, coloring
    each one with a different color as we go::

        sage: h=10; c=0.4; p=0.5;
        sage: G = Graphics()
        sage: for x in srange(1,h+1):
        ...        l = [[0,x*sqrt(3)],[-x/2,-x*sqrt(3)/2],[x/2,-x*sqrt(3)/2],[0,x*sqrt(3)]]
        ...        G+=line(l,color=hue(c + p*(x/h)))
        sage: G.show(figsize=[5,5])

    TESTS:

    From trac #4604, ensure Graphics can handle 3d objects::

        sage: g = Graphics()
        sage: g += sphere((1, 1, 1), 2)
        sage: g.show()
    """

    def __init__(self):
        """
        Create a new empty Graphics objects with all the defaults.

        EXAMPLES::

            sage: G = Graphics()
        """
        self.__fontsize = 10
        self.__show_axes = True
        self.__show_legend = False
        self.__legend_opts = {}
        self.__axes_color = (0, 0, 0)
        self.__axes_label_color = (0, 0, 0)
        self.__tick_label_color = (0, 0, 0)
        self.__axes_width = 0.8
        self.__objects = []
        self._extra_kwds = {}
        self.__bbox_extra_artists = []

    def set_aspect_ratio(self, ratio):
        """
        Set the aspect ratio, which is the ratio of height and width
        of a unit square (i.e., height/width of a unit square), or
        'automatic' (expand to fill the figure).

        INPUT:


        -  ``ratio`` - a positive real number or 'automatic'


        EXAMPLES: We create a plot of the upper half of a circle, but it
        doesn't look round because the aspect ratio is off::

            sage: P = plot(sqrt(1-x^2),(x,-1,1)); P

        So we set the aspect ratio and now it is round::

            sage: P.set_aspect_ratio(1)
            sage: P.aspect_ratio()
            1.0
            sage: P

        Note that the aspect ratio is inherited upon addition (which takes
        the max of aspect ratios of objects whose aspect ratio has been
        set)::

            sage: P + plot(sqrt(4-x^2),(x,-2,2))

        In the following example, both plots produce a circle that looks
        twice as tall as wide::

            sage: Q = circle((0,0), 0.5); Q.set_aspect_ratio(2)
            sage: (P + Q).aspect_ratio(); P+Q
            2.0
            sage: (Q + P).aspect_ratio(); Q+P
            2.0
        """
        if ratio != 'auto' and ratio != 'automatic':
            ratio = float(ratio)
            if ratio <= 0:
                raise ValueError, "the aspect ratio must be positive or 'automatic'"
        else:
            ratio = 'automatic'
        self._extra_kwds['aspect_ratio'] = ratio

    def aspect_ratio(self):
        """
        Get the current aspect ratio, which is the ratio of height to
        width of a unit square, or 'automatic'.

        OUTPUT: a positive float (height/width of a unit square), or 'automatic'
        (expand to fill the figure).

        EXAMPLES:

        The default aspect ratio for a new blank Graphics object is 'automatic'::

            sage: P = Graphics()
            sage: P.aspect_ratio()
            'automatic'

        The aspect ratio can be explicitly set different than the object's default::

            sage: P = circle((1,1), 1)
            sage: P.aspect_ratio()
            1.0
            sage: P.set_aspect_ratio(2)
            sage: P.aspect_ratio()
            2.0
            sage: P.set_aspect_ratio('automatic')
            sage: P.aspect_ratio()
            'automatic'
        """
        return self._extra_kwds.get('aspect_ratio', 'automatic')

    def legend(self, show=None):
        r"""
        Set whether or not the legend is shown by default.

        INPUT:

        -  ``show`` - (default: None) a boolean

        If called with no input, return the current legend setting.

        EXAMPLES:

        By default no legend is displayed::

            sage: P = plot(sin)
            sage: P.legend()
            False

        But if we put a label then the legend is shown::

            sage: P = plot(sin, legend_label='sin')
            sage: P.legend()
            True

        We can turn it on or off::

            sage: P.legend(False)
            sage: P.legend()
            False
            sage: P.legend(True)
            sage: P # show with the legend
        """
        if show is None:
            return self.__show_legend
        else:
            self.__show_legend = bool(show)

    def set_legend_options(self, **kwds):
        r"""
        Set various legend options.

        INPUT:

        - ``title`` - (default: None) string, the legend title

        - ``ncol`` - (default: 1) positive integer, the number of columns

        - ``columnspacing`` - (default: None) the spacing between columns

        - ``borderaxespad`` - (default: None) float, length between the axes and the legend

        - ``back_color`` - (default: (0.9, 0.9, 0.9)) This parameter can be a string
          denoting a color or an RGB tuple. The string can be a color name
          as in ('red', 'green', 'yellow', ...) or a floating point number
          like '0.8' which gets expanded to (0.8, 0.8, 0.8). The
          tuple form is just a floating point RGB tuple with all values ranging
          from 0 to 1.

        - ``handlelength`` - (default: 0.05) float, the length of the legend handles

        - ``handletextpad`` - (default: 0.5) float, the pad between the legend handle and text

        - ``labelspacing`` - (default: 0.02) float, vertical space between legend entries

        - ``loc`` - (default: 'best') May be a string, an integer or a tuple. String or
              integer inputs must be one of the following:

          - 0, 'best'

          - 1, 'upper right'

          - 2, 'upper left'

          - 3, 'lower left'

          - 4, 'lower right'

          - 5, 'right'

          - 6, 'center left'

          - 7, 'center right'

          - 8, 'lower center'

          - 9, 'upper center'

          - 10, 'center'

          - Tuple arguments represent an absolute (x, y) position on the plot
            in axes coordinates (meaning from 0 to 1 in each direction).

        - ``markerscale`` - (default: 0.6) float, how much to scale the markers in the legend.

        - ``numpoints`` - (default: 2) integer, the number of points in the legend for line

        - ``borderpad`` - (default: 0.6) float, the fractional whitespace inside the legend border
          (between 0 and 1)

        - ``font_family`` - (default: 'sans-serif') string, one of 'serif', 'sans-serif',
          'cursive', 'fantasy', 'monospace'

        - ``font_style`` - (default: 'normal') string, one of 'normal', 'italic', 'oblique'

        - ``font_variant`` - (default: 'normal') string, one of 'normal', 'small-caps'

        - ``font_weight`` - (default: 'medium') string, one of 'black', 'extra bold', 'bold',
          'semibold', 'medium', 'normal', 'light'

        - ``font_size`` - (default: 'medium') string, one of 'xx-small', 'x-small', 'small',
          'medium', 'large', 'x-large', 'xx-large' or an absolute font size (e.g. 12)

        -  ``shadow`` - (default: False) boolean - draw a shadow behind the legend

        - ``fancybox`` - (default: False) a boolean.  If True, draws a frame with a round
          fancybox.

        These are all keyword arguments.

        OUTPUT: a dictionary of all current legend options

        EXAMPLES:

        By default, no options are set::

            sage: p = plot(tan, legend_label='tan')
            sage: p.set_legend_options()
            {}

        We build a legend with a shadow::

            sage: p.set_legend_options(shadow=True)
            sage: p.set_legend_options()['shadow']
            True

        To set the legend position to the center of the plot, all these
        methods are roughly equivalent::

            sage: p.set_legend_options(loc='center'); p

        ::

            sage: p.set_legend_options(loc=10); p

        ::

            sage: p.set_legend_options(loc=(0.5,0.5)); p # aligns the bottom of the box to the center
        """
        if len(kwds) == 0:
            return self.__legend_opts
        else:
            self.__legend_opts.update(kwds)


    def get_axes_range(self):
        """
        Returns a dictionary of the range of the axes for this graphics
        object.  This is fall back to the ranges in get_minmax_data() for
        any value which the user has not explicitly set.

        .. warning::

           Changing the dictionary returned by this function does not
           change the axes range for this object.  To do that, use the
           :meth:`set_axes_range` method.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: list(sorted(L.get_axes_range().items()))
            [('xmax', 3.0), ('xmin', 1.0), ('ymax', 5.0), ('ymin', -4.0)]
            sage: L.set_axes_range(xmin=-1)
            sage: list(sorted(L.get_axes_range().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 5.0), ('ymin', -4.0)]
        """
        axes_range = self.get_minmax_data()
        axes_range.update(self._get_axes_range_dict())
        return axes_range

    def set_axes_range(self, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Set the ranges of the `x` and `y` axes.

        INPUT:


        -  ``xmin, xmax, ymin, ymax`` - floats


        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.set_axes_range(-1, 20, 0, 2)
            sage: d = L.get_axes_range()
            sage: d['xmin'], d['xmax'], d['ymin'], d['ymax']
            (-1.0, 20.0, 0.0, 2.0)
        """
        l = locals()
        axes_range = self._get_axes_range_dict()
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            if l[name] is not None:
                axes_range[name] = float(l[name])

    axes_range = set_axes_range

    def _get_axes_range_dict(self):
        """
        Returns the underlying dictionary used to store the user's
        custom ranges for the axes on this object.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L._get_axes_range_dict()
            {}
            sage: L.set_axes_range(xmin=-1)
            sage: L._get_axes_range_dict()
            {'xmin': -1.0}
        """
        try:
            return self.__axes_range
        except AttributeError:
            self.__axes_range = {}
            return self.__axes_range

    def fontsize(self, s=None):
        """
        Set the font size of axes labels and tick marks.

        INPUT:


        -  ``s`` - integer, a font size in points.


        If called with no input, return the current fontsize.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.fontsize()
            10
            sage: L.fontsize(20)
            sage: L.fontsize()
            20

        All the numbers on the axes will be very large in this plot::

            sage: L
        """
        if s is None:
            try:
                return self.__fontsize
            except AttributeError:
                self.__fontsize = 10
                return self.__fontsize
        self.__fontsize = int(s)

    def axes(self, show=None):
        """
        Set whether or not the `x` and `y` axes are shown
        by default.

        INPUT:


        -  ``show`` - bool


        If called with no input, return the current axes setting.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])

        By default the axes are displayed.

        ::

            sage: L.axes()
            True

        But we turn them off, and verify that they are off

        ::

            sage: L.axes(False)
            sage: L.axes()
            False

        Displaying L now shows a triangle but no axes.

        ::

            sage: L
        """
        if show is None:
            try:
                return self.__show_axes
            except AttributeError:
                self.__show_axes = True
                return self.__show_axes
        self.__show_axes = bool(show)

    def axes_color(self, c=None):
        """
        Set the axes color.

        If called with no input, return the current axes_color setting.

        INPUT:


        -  ``c`` - an RGB color 3-tuple, where each tuple entry
           is a float between 0 and 1


        EXAMPLES: We create a line, which has like everything a default
        axes color of black.

        ::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.axes_color()
            (0, 0, 0)

        We change the axes color to red and verify the change.

        ::

            sage: L.axes_color((1,0,0))
            sage: L.axes_color()
            (1.0, 0.0, 0.0)

        When we display the plot, we'll see a blue triangle and bright red
        axes.

        ::

            sage: L
        """
        if c is None:
            try:
                return self.__axes_color

            except AttributeError:
                self.__axes_color = (0.0, 0.0, 0.0)
                return self.__axes_color
        self.__axes_color = rgbcolor(c)

    def axes_labels(self, l=None):
        """
        Set the axes labels.

        INPUT:


        -  ``l`` - (default: None) a list of two strings or
           None


        OUTPUT: a 2-tuple of strings

        If l is None, returns the current ``axes_labels``,
        which is itself by default None. The default labels are both
        empty.

        EXAMPLES: We create a plot and put x and y axes labels on it.

        ::

            sage: p = plot(sin(x), (x, 0, 10))
            sage: p.axes_labels(['$x$','$y$'])
            sage: p.axes_labels()
            ('$x$', '$y$')

        Now when you plot p, you see x and y axes labels::

            sage: p

        Notice that some may prefer axes labels which are not
        typeset::

            sage: plot(sin(x), (x, 0, 10), axes_labels=['x','y'])
        """
        if l is None:
            try:
                return self.__axes_labels
            except AttributeError:
                self.__axes_labels = None
                return self.__axes_labels
        if not isinstance(l, (list, tuple)):
            raise TypeError, "l must be a list or tuple"
        if len(l) != 2:
            raise ValueError, "l must have length 2"
        self.__axes_labels = (str(l[0]), str(l[1]))

    def axes_label_color(self, c=None):
        r"""
        Set the color of the axes labels.

        The axes labels are placed at the edge of the x and y axes, and are
        not on by default (use the ``axes_labels`` command to
        set them; see the example below). This function just changes their
        color.

        INPUT:


        -  ``c`` - an RGB 3-tuple of numbers between 0 and 1


        If called with no input, return the current axes_label_color
        setting.

        EXAMPLES: We create a plot, which by default has axes label color
        black.

        ::

            sage: p = plot(sin, (-1,1))
            sage: p.axes_label_color()
            (0, 0, 0)

        We change the labels to be red, and confirm this::

            sage: p.axes_label_color((1,0,0))
            sage: p.axes_label_color()
            (1.0, 0.0, 0.0)

        We set labels, since otherwise we won't see anything.

        ::

            sage: p.axes_labels(['$x$ axis', '$y$ axis'])

        In the plot below, notice that the labels are red::

            sage: p
        """
        if c is None:
            try:
                return self.__axes_label_color
            except AttributeError:
                self.__axes_label_color = (0, 0, 0)
                return self.__axes_label_color
        self.__axes_label_color = rgbcolor(c)


    def axes_width(self, w=None):
        r"""
        Set the axes width. Use this to draw a plot with really fat or
        really thin axes.

        INPUT:


        -  ``w`` - a float


        If called with no input, return the current
        ``axes_width`` setting.

        EXAMPLE: We create a plot, see the default axes width (with funny
        Python float rounding), then reset the width to 10 (very fat).

        ::

            sage: p = plot(cos, (-3,3))
            sage: p.axes_width()
            0.8
            sage: p.axes_width(10)
            sage: p.axes_width()
            10.0

        Finally we plot the result, which is a graph with very fat axes.

        ::

            sage: p
        """
        if w is None:
            try:
                return self.__axes_width
            except AttributeError:
                self.__axes_width = True
                return self.__axes_width
        self.__axes_width = float(w)

    def tick_label_color(self, c=None):
        """
        Set the color of the axes tick labels.

        INPUT:


        -  ``c`` - an RGB 3-tuple of numbers between 0 and 1


        If called with no input, return the current tick_label_color
        setting.

        EXAMPLES::

            sage: p = plot(cos, (-3,3))
            sage: p.tick_label_color()
            (0, 0, 0)
            sage: p.tick_label_color((1,0,0))
            sage: p.tick_label_color()
            (1.0, 0.0, 0.0)
            sage: p
        """
        if c is None:
            try:
                return self.__tick_label_color
            except AttributeError:
                self.__tick_label_color = (0, 0, 0)
                return self.__tick_label_color
        self.__tick_label_color = rgbcolor(c)

    def _repr_(self):
        r"""
        Show this graphics objects.

        If the ``show_default`` function has been called with
        True (the default), then you'll see this graphics object displayed.
        Otherwise you'll see a text representation of it.

        EXAMPLES: We create a plot and call ``_repr_`` on it,
        which causes it to be displayed as a plot::

            sage: P = plot(cos, (-1,1))
            sage: P._repr_()
            ''

        Just doing this also displays the plot::

            sage: P

        Note that printing P with the ``print`` statement does
        not display the plot::

            sage: print P
            Graphics object consisting of 1 graphics primitive

        Now we turn off showing plots by default::

            sage: show_default(False)

        Now we just get a string. To show P you would have to do
        ``show(P)``.

        ::

            sage: P._repr_()
            'Graphics object consisting of 1 graphics primitive'
            sage: P
            Graphics object consisting of 1 graphics primitive

        Finally, we turn ``show_default`` back on::

            sage: show_default(True)
        """
        if SHOW_DEFAULT:
            self.show()
            return ''
        else:
            return self.__str__()

    def __str__(self):
        r"""
        Return string representation of this plot.

        EXAMPLES::

            sage: S = circle((0,0), 2); S.__str__()
            'Graphics object consisting of 1 graphics primitive'
            sage: print S
            Graphics object consisting of 1 graphics primitive

        .. warning::

           ``__str__`` is not called when printing lists of graphics
           objects, which can be confusing, since they will all pop
           up. One workaround is to call ``show_default``:

        For example, below when we do ``print v`` two plots are
        displayed::

            sage: v = [circle((0,0), 2), circle((2,3), 1)]
            sage: print v
            [, ]

        However, if we call ``show_default`` then we see the
        text representations of the graphics::

            sage: show_default(False)
            sage: print v
            [Graphics object consisting of 1 graphics primitive, Graphics object consisting of 1 graphics primitive]
            sage: v
            [Graphics object consisting of 1 graphics primitive,
             Graphics object consisting of 1 graphics primitive]

        ::

            sage: show_default(True)
        """
        pr, i = '', 0
        for x in self:
            pr += '\n\t%s -- %s'%(i, x)
            i += 1
        s = "Graphics object consisting of %s graphics primitives"%(len(self))
        if len(self) == 1:
            s = s[:-1]
        return s

    def __getitem__(self, i):
        """
        Returns the ith graphics primitive object:

        EXAMPLE::

            sage: G = circle((1,1),2) + circle((2,2),5); print G
            Graphics object consisting of 2 graphics primitives
            sage: G[1]
            Circle defined by (2.0,2.0) with r=5.0
        """
        return self.__objects[i]

    def __len__(self):
        """
        If G is of type Graphics, then len(G) gives the number of distinct
        graphics primitives making up that object.

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
        """
        return len(self.__objects)

    def __delitem__(self, i):
        """
        If G is of type Graphics, then del(G[i]) removes the ith distinct
        graphic primitive making up that object.

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
            sage: del(G[2])
            sage: print G
            Graphics object consisting of 2 graphics primitives
            sage: len(G)
            2
        """
        del self.__objects[int(i)]

    def __setitem__(self, i, x):
        """
        You can replace a GraphicPrimitive (point, line, circle, etc...) in
        a Graphics object G with any other GraphicPrimitive

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives

        ::

            sage: p = polygon([[1,3],[2,-2],[1,1],[1,3]]); print p
            Graphics object consisting of 1 graphics primitive

        ::

            sage: G[1] = p[0]
            sage: G    # show the plot
        """
        from sage.plot.primitive import GraphicPrimitive
        if not isinstance(x, GraphicPrimitive):
            raise TypeError, "x must be a GraphicPrimitive"
        self.__objects[int(i)] = x

    def __radd__(self, other):
        """
        Compute and return other + this graphics object.

        This only works when other is a Python int equal to 0. In all other
        cases a TypeError is raised. The main reason for this function is
        to make summing a list of graphics objects easier.

        EXAMPLES::

            sage: S = circle((0,0), 2)
            sage: print int(0) + S
            Graphics object consisting of 1 graphics primitive
            sage: print S + int(0)
            Graphics object consisting of 1 graphics primitive

        The following would fail were it not for this function::

            sage: v = [circle((0,0), 2), circle((2,3), 1)]
            sage: print sum(v)
            Graphics object consisting of 2 graphics primitives
        """
        if isinstance(other, (int, long)) and other == 0:
            return self
        raise TypeError

    def __add__(self, other):
        """
        If you have any Graphics object G1, you can always add any other
        amount of Graphics objects G2,G3,... to form a new Graphics object:
        G4 = G1 + G2 + G3.

        The xmin, xmax, ymin, and ymax properties of the graphics objects
        are expanded to include all objects in both scenes. If the aspect
        ratio property of either or both objects are set, then the larger
        aspect ratio is chosen, with 'automatic' being overridden by a
        numeric aspect ratio.

        If one of the graphics object is set to show a legend, then the
        resulting object will also be set to show a legend.  None of the
        legend options are carried over.

        EXAMPLES::

            sage: g1 = plot(abs(sqrt(x^3-1)), (x,1,5), frame=True)
            sage: g2 = plot(-abs(sqrt(x^3-1)), (x,1,5), color='red')
            sage: g1 + g2  # displays the plot

        TESTS:

        Extra keywords to show are propagated::

            sage: (g1 + g2)._extra_kwds=={'aspect_ratio': 'automatic', 'frame': True}
            True
            sage: g1.set_aspect_ratio(2)
            sage: (g1+g2).aspect_ratio()
            2.0
            sage: g2.set_aspect_ratio(3)
            sage: (g1+g2).aspect_ratio()
            3.0
        """
        if isinstance(other, int) and other == 0:
            return self
        if not isinstance(other, Graphics):
            from sage.plot.plot3d.base import Graphics3d
            if isinstance(other, Graphics3d):
                return self.plot3d() + other
            raise TypeError, "other (=%s) must be a Graphics objects"%other
        g = Graphics()
        g.__objects = self.__objects + other.__objects
        g.__show_legend = self.__show_legend or other.__show_legend
        g._extra_kwds.update(self._extra_kwds)
        g._extra_kwds.update(other._extra_kwds)
        if self.aspect_ratio()=='automatic':
            g.set_aspect_ratio(other.aspect_ratio())
        elif other.aspect_ratio()=='automatic':
            g.set_aspect_ratio(self.aspect_ratio())
        else:
            g.set_aspect_ratio( max(self.aspect_ratio(), other.aspect_ratio()))
        return g

    def add_primitive(self, primitive):
        """
        Adds a primitive to this graphics object.
        """
        self.__objects.append(primitive)

    def plot(self, *args, **kwds):
        """
        Draw a 2D plot of this graphics object, which just returns this
        object since this is already a 2D graphics object.

        EXAMPLES::

            sage: S = circle((0,0), 2)
            sage: S.plot() is S
            True
        """
        return self

    def plot3d(self, z=0, **kwds):
        """
        Returns an embedding of this 2D plot into the xy-plane of 3D space,
        as a 3D plot object. An optional parameter z can be given to
        specify the z-coordinate.

        EXAMPLES::

            sage: sum([plot(z*sin(x), 0, 10).plot3d(z) for z in range(6)]) # long time
        """
        from sage.plot.plot3d.base import Graphics3dGroup
        g = Graphics3dGroup([g.plot3d(**kwds) for g in self.__objects])
        if z:
            g = g.translate(0,0,z)
        return g

    @classmethod
    def _extract_kwds_for_show(cls, kwds, ignore=[]):
        """
        Extract keywords relevant to show() from the provided dictionary.

        EXAMPLES::

            sage: kwds = {'f': lambda x: x, 'xmin': 0, 'figsize': [1,1], 'plot_points': (40, 40)}
            sage: G_kwds = Graphics._extract_kwds_for_show(kwds, ignore='xmin')
            sage: kwds # Note how this action modifies the passed dictionary
            {'xmin': 0, 'plot_points': (40, 40), 'f': <function <lambda> at ...>}
            sage: G_kwds
            {'figsize': [1, 1]}

        This method is intended to be used with _set_extra_kwds(). Here is an
        idiom to ensure the correct keywords will get passed on to show()::

            sage: options = {} # Usually this will come from an argument
            sage: g = Graphics()
            sage: g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
        """
        result = {}
        for option in cls.SHOW_OPTIONS:
            if option not in ignore:
                try:
                    result[option] = kwds.pop(option)
                except KeyError:
                    pass
        return result

    def _set_extra_kwds(self, kwds):
        """
        Set a dictionary of keywords that will get passed on to show().

        TESTS::

            sage: g = Graphics()
            sage: g._extra_kwds
            {}
            sage: g._set_extra_kwds({'figsize': [10,10]})
            sage: g._extra_kwds
            {'figsize': [10, 10]}
            sage: g.show() # Now the (blank) plot will be extra large
        """
        self._extra_kwds = kwds

    # This dictionary has the default values for the keywords to show(). When
    # show is invoked with keyword arguments, those arguments are merged with
    # this dictionary to create a set of keywords with the defaults filled in.
    # Then, those keywords are passed on to save().

    # NOTE: If you intend to use a new parameter in show(), you should update
    # this dictionary to contain the default value for that parameter.

    SHOW_OPTIONS = dict(xmin=None, xmax=None, ymin=None, ymax=None,
                        figsize=None, fig_tight=True,
                        filename=None,
                        dpi=DEFAULT_DPI, axes=None, axes_labels=None,frame=False,
                        fontsize=None,
                        aspect_ratio=None,
                        gridlines=None, gridlinesstyle=None,
                        vgridlinesstyle=None, hgridlinesstyle=None,transparent=False,
                        show_legend=None, legend_options={},
                        axes_pad=.02, ticks_integer=False,
                        ticks=None, tick_formatter=None)

    @suboptions('legend', numpoints=2, borderpad=0.6, markerscale=0.6, shadow=False,
                labelspacing=0.02, handlelength=0.05, handletextpad=0.5, borderaxespad=None,
                loc='best', font_size='medium', font_family='sans-serif', font_style='normal',
                font_weight='medium', font_variant='normal', back_color=(0.9, 0.9, 0.9),
                title=None, ncol=1, columnspacing=None, fancybox=False)
    def show(self, **kwds):
        """
        Show this graphics image with the default image viewer.

        OPTIONAL INPUT:

        - ``filename`` - (default: None) string

        - ``dpi`` - dots per inch

        - ``figsize`` - [width, height]

        - ``fig_tight`` - (default: True) whether to clip the drawing
          tightly around drawn objects.  If True, then the resulting
          image will usually not have dimensions corresponding to
          ``figsize``.  If False, the resulting image will have
          dimensions corresponding to ``figsize``.

        - ``aspect_ratio`` - the perceived height divided by the
          perceived width. For example, if the aspect ratio is set to ``1``, circles
          will look round and a unit square will appear to have sides
          of equal length, and if the aspect ratio is set ``2``, vertical units will be
          twice as long as horizontal units, so a unit square will be twice as
          high as it is wide.  If set to ``'automatic'``, the aspect ratio
          is determined by ``figsize`` and the picture fills the figure.

        - ``axes`` - (default: True)

        - ``axes_labels`` - (default: None) list (or tuple) of two
          strings; the first is used as the label for the horizontal
          axis, and the second for the vertical axis.

        - ``fontsize`` - (default: current setting -- 10) positive
          integer; used for axes labels; if you make this very large,
          you may have to increase figsize to see all labels.

        - ``frame`` - (default: False) draw a frame around the image

        - ``gridlines`` - (default: None) can be any of the following:

          - None, False: do not add grid lines.

          - True, "automatic", "major": add grid lines at major ticks of the axes.

          - "minor": add grid at major and minor ticks.

          - [xlist,ylist]: a tuple or list containing
            two elements, where xlist (or ylist) can be
            any of the following.


            - None, False: don't add horizontal (or vertical) lines.

            - True, "automatic", "major": add horizontal (or vertical) grid lines at
              the major ticks of the axes.

            - "minor": add horizontal (or vertical) grid lines at major and minor ticks of
              axes.

            - an iterable yielding numbers n or pairs (n,opts), where n
              is the coordinate of the line and opt is a dictionary of
              MATPLOTLIB options for rendering the line.


        - ``gridlinesstyle, hgridlinesstyle, vgridlinesstyle`` -
          (default: None) a dictionary of MATPLOTLIB options for the
          rendering of the grid lines, the horizontal grid lines or the
          vertical grid lines, respectively.

        - ``linkmode`` - (default: False) If True a string containing a link
            to the produced file is returned.

        - ``transparent`` - (default: False) If True, make the background transparent.

        - ``axes_pad`` - (default: 0.02) The percentage of the axis
          range that is added to each end of each axis.  This helps
          avoid problems like clipping lines because of line-width,
          etc.  To get axes that are exactly the specified limits, set
          ``axes_pad`` to zero.

        - ``ticks_integer`` - (default: False) guarantee that the ticks
          are integers (the ``ticks`` option, if specified, will
          override this)

        - ``ticks`` - A matplotlib locator for the major ticks, or
          a number. There are several options.  For more information about
          locators, type ``from matplotlib import ticker`` and then
          ``ticker?``.

          - If this is a locator object, then it is the locator for
            the horizontal axis.  A value of None means use the default
            locator.

          - If it is a list of two locators, then the first is for the
            horizontal axis and one for the vertical axis.  A value of
            None means use the default locator (so a value of
            [None, my_locator] uses my_locator for the vertical axis and
            the default for the horizontal axis).

          - If in either case above one of the entries is a number `m`
            (something which can be coerced to a float), it will be
            replaced by a MultipleLocator which places major ticks at
            integer multiples of `m`.  See examples.

          - If in either case above one of the entries is a list of
            numbers, it will be replaced by a FixedLocator which places
            ticks at the locations specified.  This includes the case of
            of the empty list, which will give no ticks.  See examples.

        - ``tick_formatter`` - A matplotlib formatter for the major
          ticks. There are several options.  For more information about
          formatters, type ``from matplotlib import ticker`` and then
          ``ticker?``.

          If the value of this keyword is a single item, then this will
          give the formatting for the horizontal axis *only* (except for
          the ``"latex"`` option).  If it is a list or tuple, the first
          is for the horizontal axis, the second for the vertical axis.
          The options are below:

          - If one of the entries is a formatter object, then it used.
            A value of None means to use the default locator (so using
            ``tick_formatter=[None, my_formatter]`` uses my_formatter
            for the vertical axis and the default for the horizontal axis).

          - If one of the entries is a symbolic constant such as `\pi`,
            `e`, or `sqrt(2)`, ticks will be formatted nicely at rational
            multiples of this constant.

          .. warning:: This should only be used with the ``ticks`` option
             using nice rational multiples of that constant!

          - If one of the entries is the string ``"latex"``, then the
            formatting will be nice typesetting of the ticks.  This is
            intended to be used when the tick locator for at least one of
            the axes is a list including some symbolic elements.  See examples.

        - ``show_legend`` - (default: None) If True, show the legend

        - ``legend_*`` - all the options valid for :meth:`set_legend_options` prefixed with ``legend_``

        EXAMPLES::

            sage: c = circle((1,1), 1, color='red')
            sage: c.show(xmin=-1, xmax=3, ymin=-1, ymax=3)

        You could also just make the picture larger by changing ``figsize``::

            sage: c.show(figsize=8, xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can turn off the drawing of the axes::

            sage: show(plot(sin,-4,4), axes=False)

        You can also label the axes.  Putting something in dollar
        signs formats it as a mathematical expression::

            sage: show(plot(sin,-4,4), axes_labels=('$x$','$y$'))

        You can turn on the drawing of a frame around the plots::

            sage: show(plot(sin,-4,4), frame=True)

        You can make the background transparent::

            sage: plot(sin(x), (x, -4, 4), transparent=True)

        Add grid lines at the major ticks of the axes.

        ::

            sage: c = circle((0,0), 1)
            sage: c.show(gridlines=True)
            sage: c.show(gridlines="automatic")
            sage: c.show(gridlines="major")

        Add grid lines at the major and minor ticks of the axes.

        ::

            sage: u,v = var('u v')
            sage: f = exp(-(u^2+v^2))
            sage: p = plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2))
            sage: p.show(gridlines="minor")

        Add only horizontal or vertical grid lines.

        ::

            sage: p = plot(sin,-10,20)
            sage: p.show(gridlines=[None, "automatic"])
            sage: p.show(gridlines=["minor", False])

        Add grid lines at specific positions (using lists/tuples).

        ::

            sage: x, y = var('x, y')
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3)-4*(x^2+y^2-2*x)^2, \
            ...             (x,-2,2), (y,-2,2), plot_points=1000)
            sage: p.show(gridlines=[[1,0],[-1,0,1]])

        Add grid lines at specific positions (using iterators).

        ::

            sage: def maple_leaf(t):
            ...     return (100/(100+(t-pi/2)^8))*(2-sin(7*t)-cos(30*t)/2)
            sage: p = polar_plot(maple_leaf, -pi/4, 3*pi/2, color="red",plot_points=1000) # long time
            sage: p.show(gridlines=( [-3,-2.75,..,3], xrange(-1,5,2) )) # long time

        Add grid lines at specific positions (using functions).

        ::

            sage: y = x^5 + 4*x^4 - 10*x^3 - 40*x^2 + 9*x + 36
            sage: p = plot(y, -4.1, 1.1)
            sage: xlines = lambda a,b: [z for z,m in y.roots()]
            sage: p.show(gridlines=[xlines, [0]], frame=True, axes=False)

        Change the style of all the grid lines.

        ::

            sage: b = bar_chart([-3,5,-6,11], color='red')
            sage: b.show(gridlines=([-1,-0.5,..,4],True),
            ...     gridlinesstyle=dict(color="blue", linestyle=":"))

        Change the style of the horizontal or vertical grid lines
        separately.

        ::

            sage: p = polar_plot(2 + 2*cos(x), 0, 2*pi, color=hue(0.3))
            sage: p.show(gridlines=True,
            ...     hgridlinesstyle=dict(color="orange", linewidth=1.0),
            ...     vgridlinesstyle=dict(color="blue", linestyle=":"))

        Change the style of each grid line individually.

        ::

            sage: x, y = var('x, y')
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3)-4*(x^2+y^2-2*x)^2,
            ...             (x,-2,2), (y,-2,2), plot_points=1000)
            sage: p.show(gridlines=(
            ...    [
            ...     (1,{"color":"red","linestyle":":"}),
            ...     (0,{"color":"blue","linestyle":"--"})
            ...    ],
            ...    [
            ...     (-1,{"color":"red","linestyle":":"}),
            ...     (0,{"color":"blue","linestyle":"--"}),
            ...     (1,{"color":"red","linestyle":":"}),
            ...    ]
            ...    ),
            ...    gridlinesstyle=dict(marker='x',color="black"))

        Grid lines can be added to contour plots.

        ::

            sage: f = sin(x^2 + y^2)*cos(x)*sin(y)
            sage: c = contour_plot(f, (x, -4, 4), (y, -4, 4), plot_points=100)
            sage: c.show(gridlines=True, gridlinesstyle={'linestyle':':','linewidth':1, 'color':'red'})

        Grid lines can be added to matrix plots.

        ::

            sage: M = MatrixSpace(QQ,10).random_element()
            sage: matrix_plot(M).show(gridlines=True)

        By default, Sage increases the horizontal and vertical axes
        limits by a certain percentage in all directions.  This is
        controlled by the ``axes_pad`` parameter.  Increasing the range
        of the axes helps avoid problems with lines and dots being
        clipped because the linewidth extends beyond the axes.  To get
        axes limits that are exactly what is specified, set
        ``axes_pad`` to zero.  Compare the following two examples

        ::

            sage: plot(sin(x), (x, -pi, pi),thickness=2)+point((pi, -1), pointsize=15)
            sage: plot(sin(x), (x, -pi, pi),thickness=2,axes_pad=0)+point((pi, -1), pointsize=15)

        Via matplotlib, Sage allows setting of custom ticks.  See above
        for more details.

        ::

            sage: plot(sin(pi*x), (x, -8, 8)) # Labels not so helpful
            sage: plot(sin(pi*x), (x, -8, 8), ticks=2) # Multiples of 2
            sage: plot(sin(pi*x), (x, -8, 8), ticks=[[-7,-3,0,3,7],[-1/2,0,1/2]]) # Your choices
            sage: plot(sin(pi*x), (x, -8, 8), ticks=[[],[]]) # No ticks at all!

        This can be very helpful in showing certain features of plots. ::

            sage: plot(1.5/(1+e^(-x)), (x, -10, 10)) # doesn't quite show value of inflection point

        ::

            sage: plot(1.5/(1+e^(-x)), (x, -10, 10), ticks=[None, 1.5/4]) # It's right at f(x)=0.75!

        But be careful to leave enough room for at least two major ticks, so that
        the user can tell what the scale is.

        ::

            sage: plot(x^2,(x,1,8),ticks=6)
            Traceback (most recent call last):
            ...
            ValueError: Expand the range of the independent variable to allow two multiples of your tick locator (option `ticks`).

        We can also do custom formatting if you need it.  See above for full
        details.

        ::

            sage: plot(2*x+1,(x,0,5),ticks=[[0,1,e,pi,sqrt(20)],2],tick_formatter="latex")

        This is particularly useful when setting custom ticks in multiples
        of `\pi`.

        ::

            sage: plot(sin(x),(x,0,2*pi),ticks=pi/3,tick_formatter=pi)

        But keep in mind that you will get exactly the formatting you asked
        for if you specify both formatters.  The first syntax is recommended
        for best style in that case. ::

            sage: plot(arcsin(x),(x,-1,1),ticks=[None,pi/6],tick_formatter=["latex",pi]) # Nice-looking!

        ::

            sage: plot(arcsin(x),(x,-1,1),ticks=[None,pi/6],tick_formatter=[None,pi]) # Not so nice-looking

        """

        # This option should not be passed on to save().
        linkmode = kwds.pop('linkmode', False)

        if DOCTEST_MODE:
            kwds.pop('filename', None)
            self.save(DOCTEST_MODE_FILE, **kwds)
        elif EMBEDDED_MODE:
            kwds.setdefault('filename', sage.misc.misc.graphics_filename())
            self.save(**kwds)
            if linkmode == True:
                return "<img src='cell://%s'>" % kwds['filename']
            else:
                html("<img src='cell://%s'>" % kwds['filename'])
        else:
            kwds.setdefault('filename', sage.misc.misc.tmp_filename() + '.png')
            self.save(**kwds)
            os.system('%s %s 2>/dev/null 1>/dev/null &'
                      % (sage.misc.viewer.browser(), kwds['filename']))

    def xmin(self, xmin=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.xmin()
            -1.0
            sage: g.xmin(-3)
            sage: g.xmin()
            -3.0
        """
        if xmin is None:
            return self.get_axes_range()['xmin']
        else:
            self.set_axes_range(xmin=xmin)

    def xmax(self, xmax=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.xmax()
            3.0
            sage: g.xmax(10)
            sage: g.xmax()
            10.0
        """
        if xmax is None:
            return self.get_axes_range()['xmax']
        else:
            self.set_axes_range(xmax=xmax)

    def ymin(self, ymin=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.ymin()
            1.0
            sage: g.ymin(-3)
            sage: g.ymin()
            -3.0
        """
        if ymin is None:
            return self.get_axes_range()['ymin']
        else:
            self.set_axes_range(ymin=ymin)

    def ymax(self, ymax=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.ymax()
            2.0
            sage: g.ymax(10)
            sage: g.ymax()
            10.0
        """
        if ymax is None:
            return self.get_axes_range()['ymax']
        else:
            self.set_axes_range(ymax=ymax)


    def get_minmax_data(self):
        """
        Return a dictionary whose keys give the xmin, xmax, ymin, and ymax
        data for this graphic.

        .. warning::

           The returned dictionary is mutable, but changing it does
           not change the xmin/xmax/ymin/ymax data.  The minmax data is a function
           of the primitives which make up this Graphics object.  To change the
           range of the axes, call methods :meth:`xmin`, :meth:`xmax`,
           :meth:`ymin`, :meth:`ymax`, or :meth:`set_axes_range`.

        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: list(sorted(g.get_minmax_data().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 2.0), ('ymin', 1.0)]

        Note that changing ymax doesn't change the output of get_minmax_data::

            sage: g.ymax(10)
            sage: list(sorted(g.get_minmax_data().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 2.0), ('ymin', 1.0)]
        """
        objects = self.__objects
        if objects:
            minmax_data = [o.get_minmax_data() for o in objects]
            xmin = min(d['xmin'] for d in minmax_data)
            xmax = max(d['xmax'] for d in minmax_data)
            ymin = min(d['ymin'] for d in minmax_data)
            ymax = max(d['ymax'] for d in minmax_data)
            # check for NaN's: weird thing -- only way I know to check if a float
            # is a NaN is to check if it is not equal to itself.
            if xmin!=xmin:
                xmin=0; sage.misc.misc.verbose("xmin was NaN (setting to 0)", level=0)
            if xmax!=xmax:
                xmax=0; sage.misc.misc.verbose("xmax was NaN (setting to 0)", level=0)
            if ymin!=ymin:
                ymin=0; sage.misc.misc.verbose("ymin was NaN (setting to 0)", level=0)
            if ymax!=ymax:
                ymax=0; sage.misc.misc.verbose("ymax was NaN (setting to 0)", level=0)
        else:
            xmin = xmax = ymin = ymax = 0

        if xmin == xmax:
            xmin -= 1
            xmax += 1
        if ymin == ymax:
            ymin -= 1
            ymax += 1
        return {'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax}

    def matplotlib(self, filename=None,
                   xmin=None, xmax=None, ymin=None, ymax=None,
                   figsize=None, figure=None, sub=None,
                   axes=None, axes_labels=None, fontsize=None,
                   frame=False, verify=True,
                   aspect_ratio = None,
                   gridlines=None, gridlinesstyle=None,
                   vgridlinesstyle=None, hgridlinesstyle=None,
                   show_legend=None, legend_options={},
                   axes_pad=0.02, ticks_integer=None,
                   tick_formatter=None, ticks=None):
        r"""
        Return a matplotlib figure object representing the graphic

        EXAMPLES::

            sage: c = circle((1,1),1)
            sage: print c.matplotlib()
            Figure(640x480)

        To obtain the first matplotlib axes object inside of the
        figure, you can do something like the following.

        ::

            sage: p=plot(sin(x), (x, -2*pi, 2*pi))
            sage: figure=p.matplotlib()
            sage: axes=figure.axes[0]

        For input parameters, see the documentation for the
        :meth:`show` method (this function accepts all except the
        transparent argument).

        TESTS:

        We verify that #10291 is fixed::

          sage: p = plot(sin(x), (x, -2*pi, 2*pi))
          sage: figure = p.matplotlib()
          sage: axes_range = p.get_axes_range()
          sage: figure = p.matplotlib()
          sage: axes_range2 = p.get_axes_range()
          sage: axes_range == axes_range2
          True
        """
        if not isinstance(ticks, (list, tuple)):
            ticks = (ticks, None)

        from sage.symbolic.ring import SR
        if not isinstance(tick_formatter, (list, tuple)):  # make sure both formatters typeset or both don't
            if tick_formatter == "latex" or tick_formatter in SR:
                tick_formatter = (tick_formatter, "latex")
            else:
                tick_formatter = (tick_formatter, None)

        self.set_axes_range(xmin, xmax, ymin, ymax)
        d = self.get_axes_range()
        xmin = d['xmin']
        xmax = d['xmax']
        ymin = d['ymin']
        ymax = d['ymax']

        x_pad=(xmax-xmin)*float(axes_pad)
        y_pad=(ymax-ymin)*float(axes_pad)

        xmin-=x_pad
        xmax+=x_pad
        ymin-=y_pad
        ymax+=y_pad

        global do_verify
        do_verify = verify

        if axes is None:
            axes = self.__show_axes

        from matplotlib.figure import Figure
        from matplotlib import rcParams
        self.fontsize(fontsize)
        self.axes_labels(l=axes_labels)

        if figsize is not None and not isinstance(figsize, (list, tuple)):
            default_width, default_height=rcParams['figure.figsize']
            figsize=(figsize, default_height*figsize/default_width)

        if figure is None:
            figure=Figure(figsize=figsize)

        #the incoming subplot instance
        subplot = sub
        if not subplot:
            subplot = figure.add_subplot(111)
        if aspect_ratio is None:
            aspect_ratio=self.aspect_ratio()
        if aspect_ratio == 'automatic':
            subplot.set_aspect('auto', adjustable='box')
        else:
            subplot.set_aspect(aspect_ratio, adjustable='box')
        #add all the primitives to the subplot
        for g in self.__objects:
            g._render_on_subplot(subplot)
            if hasattr(g, '_bbox_extra_artists'):
                self.__bbox_extra_artists.extend(g._bbox_extra_artists)

        #add the legend if requested
        if show_legend is None:
            show_legend = self.__show_legend

        if show_legend:
            from matplotlib.font_manager import FontProperties
            lopts = dict()
            lopts.update(legend_options)
            lopts.update(self.__legend_opts)
            prop = FontProperties(family=lopts.pop('font_family'), weight=lopts.pop('font_weight'), \
                    size=lopts.pop('font_size'), style=lopts.pop('font_style'), variant=lopts.pop('font_variant'))
            color = lopts.pop('back_color')
            leg = subplot.legend(prop=prop, **lopts)
            if leg is None:
                sage.misc.misc.warn("legend requested but no items are labeled")
            else:
                # color
                lframe = leg.get_frame()
                lframe.set_facecolor(color)


        subplot.set_xlim([xmin, xmax])
        subplot.set_ylim([ymin,ymax])

        locator_options=dict(nbins=9,steps=[1,2,5,10],integer=ticks_integer)


        if axes is None:
            axes = self.__show_axes

        for spine in subplot.spines.values():
            spine.set_color(self.__axes_color)
            spine.set_linewidth(self.__axes_width)


        if frame:
            # For now, set the formatter to the old one, since that is
            # sort of what we are used to.  We should eventually look at
            # the default one to see if we like it better.

            from matplotlib.ticker import OldScalarFormatter, MaxNLocator, MultipleLocator, FixedLocator, NullLocator, Locator
            x_locator, y_locator = ticks
            if x_locator is None:
                x_locator = MaxNLocator(**locator_options)
            elif isinstance(x_locator,Locator):
                pass
            elif x_locator == []:
                x_locator = NullLocator()
            elif isinstance(x_locator,list):
                x_locator = FixedLocator(x_locator)
            else: # x_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(xmax/x_locator)-ceil(xmin/x_locator)>1:
                    x_locator=MultipleLocator(float(x_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the independent variable to allow two multiples of your tick locator (option `ticks`).')
            if y_locator is None:
                y_locator = MaxNLocator(**locator_options)
            elif isinstance(y_locator,Locator):
                pass
            elif y_locator == []:
                y_locator = NullLocator()
            elif isinstance(y_locator,list):
                y_locator = FixedLocator(y_locator)
            else: # y_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(ymax/y_locator)-ceil(ymin/y_locator)>1:
                    y_locator=MultipleLocator(float(y_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the dependent variable to allow two multiples of your tick locator (option `ticks`).')

            x_formatter, y_formatter = tick_formatter
            from matplotlib.ticker import FuncFormatter
            from sage.misc.latex import latex
            if x_formatter is None:
                x_formatter = OldScalarFormatter()
            elif x_formatter in SR:
                from misc import _multiple_of_constant
                x_const = x_formatter
                x_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,x_const))
            elif x_formatter == "latex":
                x_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))
            if y_formatter is None:
                y_formatter = OldScalarFormatter()
            elif y_formatter in SR:
                from misc import _multiple_of_constant
                y_const = y_formatter
                y_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,y_const))
            elif y_formatter == "latex":
                y_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))

            subplot.xaxis.set_major_locator(x_locator)
            subplot.yaxis.set_major_locator(y_locator)
            subplot.xaxis.set_major_formatter(x_formatter)
            subplot.yaxis.set_major_formatter(y_formatter)

            subplot.set_frame_on(True)
            if axes:
                if ymin<=0 and ymax>=0:
                    subplot.axhline(color=self.__axes_color,
                                    linewidth=self.__axes_width)
                if xmin<=0 and xmax>=0:
                    subplot.axvline(color=self.__axes_color,
                                    linewidth=self.__axes_width)

        elif axes:
            ymiddle=False
            xmiddle=False
            if xmin>0:
                subplot.spines['right'].set_visible(False)
                subplot.spines['left'].set_position(('outward',10))
                subplot.yaxis.set_ticks_position('left')
                subplot.yaxis.set_label_position('left')
                yaxis='left'
            elif xmax<0:
                subplot.spines['left'].set_visible(False)
                subplot.spines['right'].set_position(('outward',10))
                subplot.yaxis.set_ticks_position('right')
                subplot.yaxis.set_label_position('right')
                yaxis='right'
            else:
                subplot.spines['left'].set_position('zero')
                subplot.yaxis.set_ticks_position('left')
                subplot.yaxis.set_label_position('left')
                subplot.spines['right'].set_visible(False)
                ymiddle=True
                yaxis='left'

            if ymin>0:
                subplot.spines['top'].set_visible(False)
                subplot.spines['bottom'].set_position(('outward',10))
                subplot.xaxis.set_ticks_position('bottom')
                subplot.xaxis.set_label_position('bottom')
                xaxis='bottom'
            elif ymax<0:
                subplot.spines['bottom'].set_visible(False)
                subplot.spines['top'].set_position(('outward',10))
                subplot.xaxis.set_ticks_position('top')
                subplot.xaxis.set_label_position('top')
                xaxis='top'
            else:
                subplot.spines['bottom'].set_position('zero')
                subplot.xaxis.set_ticks_position('bottom')
                subplot.xaxis.set_label_position('bottom')
                subplot.spines['top'].set_visible(False)
                xmiddle=True
                xaxis='bottom'

            # For now, set the formatter to the old one, since that is
            # sort of what we are used to.  We should eventually look at
            # the default one to see if we like it better.

            from matplotlib.ticker import OldScalarFormatter, MaxNLocator, MultipleLocator, FixedLocator, NullLocator, Locator
            x_locator, y_locator = ticks
            if x_locator is None:
                x_locator = MaxNLocator(**locator_options)
            elif isinstance(x_locator,Locator):
                pass
            elif x_locator == []:
                x_locator = NullLocator()
            elif isinstance(x_locator,list):
                x_locator = FixedLocator(x_locator)
            else: # x_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(xmax/x_locator)-ceil(xmin/x_locator)>1:
                    x_locator=MultipleLocator(float(x_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the independent variable to allow two multiples of your tick locator (option `ticks`).')
            if y_locator is None:
                y_locator = MaxNLocator(**locator_options)
            elif isinstance(y_locator,Locator):
                pass
            elif y_locator == []:
                y_locator = NullLocator()
            elif isinstance(y_locator,list):
                y_locator = FixedLocator(y_locator)
            else: # y_locator is a number which can be made a float
                from sage.functions.other import ceil, floor
                if floor(ymax/y_locator)-ceil(ymin/y_locator)>1:
                    y_locator=MultipleLocator(float(y_locator))
                else: # not enough room for two major ticks
                    raise ValueError('Expand the range of the dependent variable to allow two multiples of your tick locator (option `ticks`).')

            x_formatter, y_formatter = tick_formatter
            from matplotlib.ticker import FuncFormatter
            from sage.misc.latex import latex
            from sage.symbolic.ring import SR
            if x_formatter is None:
                x_formatter = OldScalarFormatter()
            elif x_formatter in SR:
                from misc import _multiple_of_constant
                x_const = x_formatter
                x_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,x_const))
            elif x_formatter == "latex":
                x_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))
            if y_formatter is None:
                y_formatter = OldScalarFormatter()
            elif y_formatter in SR:
                from misc import _multiple_of_constant
                y_const = y_formatter
                y_formatter = FuncFormatter(lambda n,pos: _multiple_of_constant(n,pos,y_const))
            elif y_formatter == "latex":
                y_formatter = FuncFormatter(lambda n,pos: '$%s$'%latex(n))

            subplot.xaxis.set_major_locator(x_locator)
            subplot.yaxis.set_major_locator(y_locator)
            subplot.xaxis.set_major_formatter(x_formatter)
            subplot.yaxis.set_major_formatter(y_formatter)

            # Make ticklines go on both sides of the axes
            #             if xmiddle:
            #                 for t in subplot.xaxis.get_majorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(8)
            #                 for t in subplot.xaxis.get_minorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(4)

            #             if ymiddle:
            #                 for t in subplot.yaxis.get_majorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(8)
            #                 for t in subplot.yaxis.get_minorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(4)

            # Make the zero tick labels disappear if the axes cross
            # inside the picture
            if xmiddle and ymiddle:
                subplot.yaxis.set_major_formatter(SelectiveFormatter(subplot.yaxis.get_major_formatter(),skip_values=[0]))
                subplot.xaxis.set_major_formatter(SelectiveFormatter(subplot.xaxis.get_major_formatter(),skip_values=[0]))

        else:
            for spine in subplot.spines.values():
                spine.set_visible(False)
            from matplotlib.ticker import NullFormatter, NullLocator
            subplot.xaxis.set_major_formatter(NullFormatter())
            subplot.yaxis.set_major_formatter(NullFormatter())
            subplot.xaxis.set_major_locator(NullLocator())
            subplot.yaxis.set_major_locator(NullLocator())

        if frame or axes:
            # Make minor tickmarks, unless we specify fixed ticks or no ticks
            from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullLocator
            if isinstance(x_locator, (NullLocator, FixedLocator)):
                subplot.xaxis.set_minor_locator(NullLocator())
            else:
                subplot.xaxis.set_minor_locator(AutoMinorLocator())
            if isinstance(y_locator, (NullLocator, FixedLocator)):
                subplot.yaxis.set_minor_locator(NullLocator())
            else:
                subplot.yaxis.set_minor_locator(AutoMinorLocator())

            ticklabels=subplot.xaxis.get_majorticklabels() + \
                subplot.xaxis.get_minorticklabels() + \
                subplot.yaxis.get_majorticklabels() + \
                subplot.yaxis.get_minorticklabels()
            for ticklabel in ticklabels:
                ticklabel.set_fontsize(self.__fontsize)
                ticklabel.set_color(self.__tick_label_color)

            ticklines=subplot.xaxis.get_majorticklines() + \
                subplot.xaxis.get_minorticklines() + \
                subplot.yaxis.get_majorticklines() + \
                subplot.yaxis.get_minorticklines()
            for tickline in ticklines:
                tickline.set_color(self.__axes_color)


        if gridlines is not None:
            if isinstance(gridlines, (list, tuple)):
                vgridlines,hgridlines=gridlines
            else:
                hgridlines=gridlines
                vgridlines=gridlines

            if gridlinesstyle is None:
                # Set up the default grid style
                gridlinesstyle=dict(color='black',linestyle=':',linewidth=0.5)

            vgridstyle=gridlinesstyle.copy()
            if vgridlinesstyle is not None:
                vgridstyle.update(vgridlinesstyle)

            hgridstyle=gridlinesstyle.copy()
            if hgridlinesstyle is not None:
                hgridstyle.update(hgridlinesstyle)

            if hgridlines=='minor':
                hgridstyle['which']='both'
            if vgridlines=='minor':
                vgridstyle['which']='both'

            if hasattr(hgridlines, '__iter__'):
                hlines=iter(hgridlines)
                hgridstyle.pop("minor",None)
                for hline in hlines:
                    if isinstance(hline, (list, tuple)):
                        hl, style=hline
                        st=hgridstyle.copy()
                        st.update(style)
                    else:
                        hl=hline
                        st=hgridstyle
                    subplot.axhline(hl,**st)
            else:
                if hgridlines not in (None, False):
                    subplot.yaxis.grid(True, **hgridstyle)

            if hasattr(vgridlines, '__iter__'):
                vlines=iter(vgridlines)
                vgridstyle.pop("minor",None)
                for vline in vlines:
                    if isinstance(vline, (list, tuple)):
                        vl, style=vline
                        st=vgridstyle.copy()
                        st.update(style)
                    else:
                        vl=vline
                        st=vgridstyle
                    subplot.axvline(vl,**st)
            else:
                if vgridlines not in (None, False):
                    subplot.xaxis.grid(True, **vgridstyle)



        if self.__axes_labels is not None:
            label_options={}
            label_options['color']=self.__axes_label_color
            label_options['size']=self.__fontsize
            subplot.set_xlabel(self.__axes_labels[0], **label_options)
            subplot.set_ylabel(self.__axes_labels[1], **label_options)


            if axes is True and frame is False:
                # We set the label positions according to where we are
                # drawing the axes.
                if xaxis=='bottom':
                    yaxis_labely=subplot.get_ylim()[1]
                    yaxis_labeloffset=8
                    yaxis_vert='bottom'
                    xaxis_labely=0
                    xaxis_vert='baseline'
                else:
                    yaxis_labely=subplot.get_ylim()[0]
                    yaxis_labeloffset=-8
                    yaxis_vert='top'
                    xaxis_labely=1
                    xaxis_vert='top'

                if yaxis=='left':
                    xaxis_labelx=subplot.get_xlim()[1]
                    xaxis_labeloffset=8
                    xaxis_horiz='left'
                    yaxis_labelx=0
                else:
                    xaxis_labelx=subplot.get_xlim()[0]
                    xaxis_labeloffset=-8
                    xaxis_horiz='right'
                    yaxis_labelx=1

                from matplotlib.transforms import offset_copy
                xlabel=subplot.xaxis.get_label()
                xlabel.set_horizontalalignment(xaxis_horiz)
                xlabel.set_verticalalignment(xaxis_vert)
                trans=subplot.spines[xaxis].get_transform()
                labeltrans=offset_copy(trans, figure, x=xaxis_labeloffset, y=0, units='points')
                subplot.xaxis.set_label_coords(x=xaxis_labelx,y=xaxis_labely,transform=labeltrans)

                ylabel=subplot.yaxis.get_label()
                ylabel.set_horizontalalignment('center')
                ylabel.set_verticalalignment(yaxis_vert)
                ylabel.set_rotation('horizontal')
                trans=subplot.spines[yaxis].get_transform()
                labeltrans=offset_copy(trans, figure, x=0, y=yaxis_labeloffset, units='points')
                subplot.yaxis.set_label_coords(x=yaxis_labelx,y=yaxis_labely,transform=labeltrans)

        # This option makes the xlim and ylim limits not take effect
        # todo: figure out which limits were specified, and let the
        # free limits autoscale
        #subplot.autoscale_view(tight=True)
        return figure

    # ALLOWED_EXTENSIONS is the list of recognized formats.
    # filename argument is written explicitly so that it can be used as a
    # positional one, which is a very likely usage for this function.
    @suboptions('legend', numpoints=2, borderpad=0.6, markerscale=0.6, shadow=False,
                labelspacing=0.02, handlelength=0.05, handletextpad=0.5, borderaxespad=None,
                loc='best', font_size='medium', font_family='sans-serif', font_style='normal',
                font_weight='medium', font_variant='normal', back_color=(0.9, 0.9, 0.9),
                title=None, ncol=1, columnspacing=None, fancybox=False)
    def save(self, filename=None, **kwds):
        r"""
        Save the graphics to an image file.

        INPUT:

        - ``filename`` -- a string (default: autogenerated), the filename and
          the image format given by the extension, which can be one of the
          following:

            * ``.eps``,

            * ``.pdf``,

            * ``.png``,

            * ``.ps``,

            * ``.sobj`` (for a Sage object you can load later),

            * ``.svg``,

            * empty extension will be treated as ``.sobj``.

        All other keyword arguments will be passed to the plotter.

        OUTPUT:

        - none.

        EXAMPLES::

            sage: c = circle((1,1), 1, color='red')
            sage: filename = os.path.join(SAGE_TMP, 'test.png')
            sage: c.save(filename, xmin=-1, xmax=3, ymin=-1, ymax=3)

        To make a figure bigger or smaller, use ``figsize``::

            sage: c.save(filename, figsize=5, xmin=-1, xmax=3, ymin=-1, ymax=3)

        By default, the figure grows to include all of the graphics and text,
        so the final image may not be exactly the figure size you specified.
        If you want a figure to be exactly a certain size, specify the keyword
        ``fig_tight=False``::

            sage: c.save(filename, figsize=[8,4], fig_tight=False,
            ...       xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can also pass extra options to the plot command instead of this
        method, e.g. ::

            sage: plot(x^2 - 5, (x, 0, 5), ymin=0).save(
            ...       sage.misc.misc.tmp_filename() + '.png')

        will save the same plot as the one shown by this command::

            sage: plot(x^2 - 5, (x, 0, 5), ymin=0)

        (This test verifies that Trac #8632 is fixed.)

        TESTS:

        Legend labels should save correctly::

            sage: P = plot(x,(x,0,1),legend_label='$xyz$')
            sage: P.set_legend_options(back_color=(1,0,0))
            sage: P.set_legend_options(loc=7)
            sage: filename=os.path.join(SAGE_TMP, 'test.png')
            sage: P.save(filename)

        This plot should save with the frame shown, showing Trac #7524
        is fixed (same issue as #7981 and #8632)::

            sage: var('x,y')
            (x, y)
            sage: a = plot_vector_field((x,-y),(x,-1,1),(y,-1,1))
            sage: filename=os.path.join(SAGE_TMP, 'test2.png')
            sage: a.save(filename)
        """
        options = dict()
        options.update(self.SHOW_OPTIONS)
        options.update(self._extra_kwds)
        options.update(kwds)
        dpi = options.pop('dpi')
        transparent = options.pop('transparent')
        fig_tight = options.pop('fig_tight')

        if filename is None:
            filename = options.pop('filename')
        if filename is None:
            filename = sage.misc.misc.graphics_filename()
        ext = os.path.splitext(filename)[1].lower()

        if ext not in ALLOWED_EXTENSIONS:
            raise ValueError("allowed file extensions for images are '"
                             + "', '".join(ALLOWED_EXTENSIONS) + "'!")
        elif ext in ['', '.sobj']:
            SageObject.save(self, filename)
        else:
            figure = self.matplotlib(**options)
            # You can output in PNG, PS, EPS, PDF, or SVG format, depending on the file extension.
            # matplotlib looks at the file extension to see what the renderer should be.
            # The default is FigureCanvasAgg for PNG's because this is by far the most
            # common type of files rendered, like in the notebook, for example.
            # if the file extension is not '.png', then matplotlib will handle it.
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            figure.set_canvas(FigureCanvasAgg(figure))
            # this messes up the aspect ratio!
            #figure.canvas.mpl_connect('draw_event', pad_for_tick_labels)

            # tight_layout adjusts the *subplot* parameters so ticks aren't cut off, etc.
            figure.tight_layout()

            if fig_tight is True:
                figure.savefig(filename, dpi=dpi, bbox_inches='tight',
                    bbox_extra_artists=self.__bbox_extra_artists,
                    transparent=transparent)
            else:
                figure.savefig(filename, dpi=dpi,
                           transparent=transparent)

#Currently not used - see comment immediately above about
#figure.canvas.mpl_connect('draw_event', pad_for_tick_labels)
# TODO - figure out how to use this, add documentation
#def pad_for_tick_labels(event):
#    import matplotlib.transforms as mtransforms
#    figure=event.canvas.figure
#    bboxes = []
#    for ax in figure.axes:
#        bbox = ax.xaxis.get_label().get_window_extent()
#        # the figure transform goes from relative coords->pixels and we
#        # want the inverse of that
#        bboxi = bbox.inverse_transformed(figure.transFigure)
#        bboxes.append(bboxi)
#
#        bbox = ax.yaxis.get_label().get_window_extent()
#        bboxi = bbox.inverse_transformed(figure.transFigure)
#        bboxes.append(bboxi)
#        for label in (ax.get_xticklabels()+ax.get_yticklabels() \
#                          + ax.get_xticklabels(minor=True) \
#                          +ax.get_yticklabels(minor=True)):
#            bbox = label.get_window_extent()
#            bboxi = bbox.inverse_transformed(figure.transFigure)
#            bboxes.append(bboxi)
#
#    # this is the bbox that bounds all the bboxes, again in relative
#    # figure coords
#    bbox = mtransforms.Bbox.union(bboxes)
#    adjusted=adjust_figure_to_contain_bbox(figure,bbox)
#
#    if adjusted:
#        figure.canvas.draw()
#    return False
#
#Currently not used - see comment above about
#figure.canvas.mpl_connect('draw_event', pad_for_tick_labels)
# TODO - figure out how to use this, add documentation
#def adjust_figure_to_contain_bbox(fig, bbox,pad=1.1):
#    """
#    For each amount we are over (in axes coordinates), we adjust by over*pad
#    to give ourselves a bit of padding.
#    """
#    left=fig.subplotpars.left
#    bottom=fig.subplotpars.bottom
#    right=fig.subplotpars.right
#    top=fig.subplotpars.top
#
#    adjusted=False
#    if bbox.xmin<0:
#        left-=bbox.xmin*pad
#        adjusted=True
#    if bbox.ymin<0:
#        bottom-=bbox.ymin*pad
#        adjusted=True
#    if bbox.xmax>1:
#        right-=(bbox.xmax-1)*pad
#        adjusted=True
#    if bbox.ymax>1:
#        top-=(bbox.ymax-1)*pad
#        adjusted=True
#
#    if left<right and bottom<top:
#        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
#        return adjusted
#    else:
#        return False

_SelectiveFormatterClass = None

def SelectiveFormatter(formatter, skip_values):

    """
    This matplotlib formatter selectively omits some tick values and
    passes the rest on to a specified formatter.

    EXAMPLES:

    This example is almost straight from a matplotlib example.

    ::

        sage: from sage.plot.plot import SelectiveFormatter
        sage: import matplotlib.pyplot as plt
        sage: import numpy
        sage: fig=plt.figure()
        sage: ax=fig.add_subplot(111)
        sage: t = numpy.arange(0.0, 2.0, 0.01)
        sage: s = numpy.sin(2*numpy.pi*t)
        sage: p = ax.plot(t, s)
        sage: formatter=SelectiveFormatter(ax.xaxis.get_major_formatter(),skip_values=[0,1])
        sage: ax.xaxis.set_major_formatter(formatter)
        sage: fig.savefig(os.path.join(SAGE_TMP, 'test.png'))
    """
    global _SelectiveFormatterClass
    if _SelectiveFormatterClass is None:

        from matplotlib.ticker import Formatter

        class _SelectiveFormatterClass(Formatter):
            def __init__(self, formatter,skip_values):
                """
                Initialize a SelectiveFormatter object.

                INPUT:

                  - formatter -- the formatter object to which we should pass labels

                  - skip_values -- a list of values that we should skip when
                    formatting the tick labels

                EXAMPLES::

                    sage: from sage.plot.plot import SelectiveFormatter
                    sage: import matplotlib.pyplot as plt
                    sage: import numpy
                    sage: fig=plt.figure()
                    sage: ax=fig.add_subplot(111)
                    sage: t = numpy.arange(0.0, 2.0, 0.01)
                    sage: s = numpy.sin(2*numpy.pi*t)
                    sage: line=ax.plot(t, s)
                    sage: formatter=SelectiveFormatter(ax.xaxis.get_major_formatter(),skip_values=[0,1])
                    sage: ax.xaxis.set_major_formatter(formatter)
                    sage: fig.savefig(os.path.join(SAGE_TMP, 'test.png'))
                """
                self.formatter=formatter
                self.skip_values=skip_values
            def set_locs(self, locs):
                """
                Set the locations for the ticks that are not skipped.

                EXAMPLES::
                    sage: from sage.plot.plot import SelectiveFormatter
                    sage: import matplotlib.ticker
                    sage: formatter=SelectiveFormatter(matplotlib.ticker.Formatter(),skip_values=[0,200])
                    sage: formatter.set_locs([i*100 for i in range(10)])
                """
                self.formatter.set_locs([l for l in locs if l not in self.skip_values])
            def __call__(self, x, *args, **kwds):
                """
                Return the format for tick val *x* at position *pos*

                EXAMPLES::

                    sage: from sage.plot.plot import SelectiveFormatter
                    sage: import matplotlib.ticker
                    sage: formatter=SelectiveFormatter(matplotlib.ticker.FixedFormatter(['a','b']),skip_values=[0,2])
                    sage: [formatter(i,1) for i in range(10)]
                    ['', 'b', '', 'b', 'b', 'b', 'b', 'b', 'b', 'b']
                """
                if x in self.skip_values:
                    return ''
                else:
                    return self.formatter(x, *args, **kwds)

    return _SelectiveFormatterClass(formatter, skip_values)


def xydata_from_point_list(points):
    r"""
    Returns two lists (xdata, ydata), each coerced to a list of floats,
    which correspond to the x-coordinates and the y-coordinates of the
    points.

    The points parameter can be a list of 2-tuples or some object that
    yields a list of one or two numbers.

    This function can potentially be very slow for large point sets.

    TESTS::

        sage: from sage.plot.plot import xydata_from_point_list
        sage: xydata_from_point_list([CC(0), CC(1)])   # ticket 8082
        ([0.0, 1.0], [0.0, 0.0])

    This function should work for anything than can be turned into a
    list, such as iterators and such (see ticket #10478)::

        sage: xydata_from_point_list(iter([(0,0), (sqrt(3), 2)]))
        ([0.0, 1.7320508075688772], [0.0, 2.0])
        sage: xydata_from_point_list((x, x^2) for x in range(5))
        ([0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 4.0, 9.0, 16.0])
        sage: xydata_from_point_list(enumerate(prime_range(1, 15)))
        ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [2.0, 3.0, 5.0, 7.0, 11.0, 13.0])
        sage: from itertools import izip; xydata_from_point_list(izip([2,3,5,7], [11, 13, 17, 19]))
        ([2.0, 3.0, 5.0, 7.0], [11.0, 13.0, 17.0, 19.0])
    """
    from sage.rings.complex_number import ComplexNumber
    if not isinstance(points, (list,tuple)):
        points = list(points)
        try:
            points = [[float(z) for z in points]]
        except TypeError:
            pass
    elif len(points)==2 and not isinstance(points[0],(list,tuple,ComplexNumber)):
        try:
            points = [[float(z) for z in points]]
        except TypeError:
            pass

    if len(points)>0 and len(list(points[0]))!=2:
        raise ValueError, "points must have 2 coordinates in a 2d line"


    xdata = [float(z[0]) for z in points]
    ydata = [float(z[1]) for z in points]

    return xdata, ydata

@rename_keyword(color='rgbcolor')
@options(alpha=1, thickness=1, fill=False, fillcolor='automatic', fillalpha=0.5, rgbcolor=(0,0,1), plot_points=200,
         adaptive_tolerance=0.01, adaptive_recursion=5, detect_poles = False, exclude = None, legend_label=None,
         __original_opts=True, aspect_ratio='automatic')
def plot(funcs, *args, **kwds):
    r"""
    Use plot by writing

    ``plot(X, ...)``

    where `X` is a Sage object (or list of Sage objects) that
    either is callable and returns numbers that can be coerced to
    floats, or has a plot method that returns a
    ``GraphicPrimitive`` object.

    There are many other specialized 2D plot commands available
    in Sage, such as ``plot_slope_field``, as well as various
    graphics primitives like Arrow; type ``sage.plot.plot?`` for
    a current list.

    Type ``plot.options`` for a dictionary of the default
    options for plots. You can change this to change the defaults for
    all future plots. Use ``plot.reset()`` to reset to the
    default options.

    PLOT OPTIONS:

    - ``plot_points`` - (default: 200) the minimal number of plot points.

    - ``adaptive_recursion`` - (default: 5) how many levels of recursion to go
      before giving up when doing adaptive refinement.  Setting this to 0
      disables adaptive refinement.

    - ``adaptive_tolerance`` - (default: 0.01) how large a difference should be
      before the adaptive refinement code considers it significant.  See the
      documentation further below for more information, starting at "the
      algorithm used to insert".

    - ``xmin`` - starting x value

    - ``xmax`` - ending x value

    - ``color`` - an RGB tuple (r,g,b) with each of r,g,b between 0 and 1,
      or a color name as a string (e.g., 'purple'), or an HTML color
      such as '#aaff0b'.

    - ``detect_poles`` - (Default: False) If set to True poles are detected.
      If set to "show" vertical asymptotes are drawn.

    - ``legend_label`` - the label for this item in the legend

    APPEARANCE OPTIONS:

    The following options affect the appearance of
    the line through the points on the graph of `X` (these are
    the same as for the line function):

    INPUT:

    - ``alpha`` - How transparent the line is

    - ``thickness`` - How thick the line is

    - ``rgbcolor`` - The color as an RGB tuple

    - ``hue`` - The color given as a hue

    Any MATPLOTLIB line option may also be passed in.  E.g.,

    - ``linestyle`` - The style of the line, which is one of
       - ``"-"`` (solid) -- default
       - ``"--"`` (dashed)
       - ``"-."`` (dash dot)
       - ``":"`` (dotted)
       - ``"None"`` or ``" "`` or ``""`` (nothing)

       The linestyle can also be prefixed with a drawing style (e.g., ``"steps--"``)

       - ``"default"`` (connect the points with straight lines)
       - ``"steps"`` or ``"steps-pre"`` (step function; horizontal
         line is to the left of point)
       - ``"steps-mid"`` (step function; points are in the middle of
         horizontal lines)
       - ``"steps-post"`` (step function; horizontal line is to the
         right of point)

    - ``marker``  - The style of the markers, which is one of
       - ``"None"`` or ``" "`` or ``""`` (nothing) -- default
       - ``","`` (pixel), ``"."`` (point)
       - ``"_"`` (horizontal line), ``"|"`` (vertical line)
       - ``"o"`` (circle), ``"p"`` (pentagon), ``"s"`` (square), ``"x"`` (x), ``"+"`` (plus), ``"*"`` (star)
       - ``"D"`` (diamond), ``"d"`` (thin diamond)
       - ``"H"`` (hexagon), ``"h"`` (alternative hexagon)
       - ``"<"`` (triangle left), ``">"`` (triangle right), ``"^"`` (triangle up), ``"v"`` (triangle down)
       - ``"1"`` (tri down), ``"2"`` (tri up), ``"3"`` (tri left), ``"4"`` (tri right)
       - ``0`` (tick left), ``1`` (tick right), ``2`` (tick up), ``3`` (tick down)
       - ``4`` (caret left), ``5`` (caret right), ``6`` (caret up), ``7`` (caret down)
       - ``"$...$"`` (math TeX string)

    - ``markersize`` - the size of the marker in points

    - ``markeredgecolor`` -- the color of the marker edge

    - ``markerfacecolor`` -- the color of the marker face

    - ``markeredgewidth`` - the size of the marker edge in points

    - ``exclude`` - (Default: None) values which are excluded from the plot range.
      Either a list of real numbers, or an equation in one variable.

    FILLING OPTIONS:

    - ``fill`` - (Default: False) One of:

      - "axis" or True: Fill the area between the function and the x-axis.

      - "min": Fill the area between the function and its minimal value.

      - "max": Fill the area between the function and its maximal value.

      - a number c: Fill the area between the function and the horizontal line y = c.

      - a function g: Fill the area between the function that is plotted and g.

      - a dictionary d (only if a list of functions are plotted):
        The keys of the dictionary should be integers.
        The value of d[i] specifies the fill options for the i-th function in the list.
        If d[i] == [j]: Fill the area between the i-th and the j-th function in the list.
        (But if d[i] == j: Fill the area between the i-th function in the list and the
        horizontal line y = j.)

    - ``fillcolor`` - (default: 'automatic') The color of the fill.
      Either 'automatic' or a color.

    - ``fillalpha`` - (default: 0.5) How transparent the fill is.
      A number between 0 and 1.

    Note that this function does NOT simply sample equally spaced
    points between xmin and xmax. Instead it computes equally spaced
    points and add small perturbations to them. This reduces the
    possibility of, e.g., sampling sin only at multiples of
    `2\pi`, which would yield a very misleading graph.

    EXAMPLES: We plot the sin function::

        sage: P = plot(sin, (0,10)); print P
        Graphics object consisting of 1 graphics primitive
        sage: len(P)     # number of graphics primitives
        1
        sage: len(P[0])  # how many points were computed (random)
        225
        sage: P          # render

    ::

        sage: P = plot(sin, (0,10), plot_points=10); print P
        Graphics object consisting of 1 graphics primitive
        sage: len(P[0])  # random output
        32
        sage: P          # render

    We plot with ``randomize=False``, which makes the initial sample points
    evenly spaced (hence always the same). Adaptive plotting might
    insert other points, however, unless ``adaptive_recursion=0``.

    ::

        sage: p=plot(1, (x,0,3), plot_points=4, randomize=False, adaptive_recursion=0)
        sage: list(p[0])
        [(0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]

    Some colored functions::

        sage: plot(sin, 0, 10, color='purple')
        sage: plot(sin, 0, 10, color='#ff00ff')

    We plot several functions together by passing a list of functions
    as input::

        sage: plot([sin(n*x) for n in [1..4]], (0, pi))

    We can also build a plot step by step from an empty plot::

        sage: a = plot([]); a       # passing an empty list returns an empty plot (Graphics() object)
        sage: a += plot(x**2); a    # append another plot
        sage: a += plot(x**3); a    # append yet another plot


    The function `\sin(1/x)` wiggles wildly near `0`.
    Sage adapts to this and plots extra points near the origin.

    ::

        sage: plot(sin(1/x), (x, -1, 1))

    Via the matplotlib library, Sage makes it easy to tell whether
    a graph is on both sides of both axes, as the axes only cross
    if the origin is actually part of the viewing area::

        sage: plot(x^3,(x,0,2))  # this one has the origin
        sage: plot(x^3,(x,1,2))  # this one does not

    Another thing to be aware of with axis labeling is that when
    the labels have quite different orders of magnitude or are very
    large, scientific notation (the `e` notation for powers of ten) is used::

        sage: plot(x^2,(x,480,500))  # this one has no scientific notation
        sage: plot(x^2,(x,300,500))  # this one has scientific notation on y-axis

    You can put a legend with ``legend_label`` (the legend is only put
    once in the case of multiple functions)::

        sage: plot(exp(x), 0, 2, legend_label='$e^x$')

    Sage understands TeX, so these all are slightly different, and you can choose
    one based on your needs::

        sage: plot(sin, legend_label='sin')
        sage: plot(sin, legend_label='$sin$')
        sage: plot(sin, legend_label='$\sin$')

    Note that the independent variable may be omitted if there is no
    ambiguity::

        sage: plot(sin(1/x), (-1, 1))

    The algorithm used to insert extra points is actually pretty
    simple. On the picture drawn by the lines below::

        sage: p = plot(x^2, (-0.5, 1.4)) + line([(0,0), (1,1)], color='green')
        sage: p += line([(0.5, 0.5), (0.5, 0.5^2)], color='purple')
        sage: p += point(((0, 0), (0.5, 0.5), (0.5, 0.5^2), (1, 1)), color='red', pointsize=20)
        sage: p += text('A', (-0.05, 0.1), color='red')
        sage: p += text('B', (1.01, 1.1), color='red')
        sage: p += text('C', (0.48, 0.57), color='red')
        sage: p += text('D', (0.53, 0.18), color='red')
        sage: p.show(axes=False, xmin=-0.5, xmax=1.4, ymin=0, ymax=2)

    You have the function (in blue) and its approximation (in green)
    passing through the points A and B. The algorithm finds the
    midpoint C of AB and computes the distance between C and D. If that
    distance exceeds the ``adaptive_tolerance`` threshold (*relative* to
    the size of the initial plot subintervals), the point D is
    added to the curve.  If D is added to the curve, then the
    algorithm is applied recursively to the points A and D, and D and
    B. It is repeated ``adaptive_recursion`` times (5, by default).

    The actual sample points are slightly randomized, so the above
    plots may look slightly different each time you draw them.

    We draw the graph of an elliptic curve as the union of graphs of 2
    functions.

    ::

        sage: def h1(x): return abs(sqrt(x^3  - 1))
        sage: def h2(x): return -abs(sqrt(x^3  - 1))
        sage: P = plot([h1, h2], 1,4)
        sage: P          # show the result

    We can also directly plot the elliptic curve::

        sage: E = EllipticCurve([0,-1])
        sage: plot(E, (1, 4), color=hue(0.6))

    We can change the line style as well::

        sage: plot(sin(x), (x, 0, 10), linestyle='-.')

    If we have an empty linestyle and specify a marker, we can see the
    points that are actually being plotted::

        sage: plot(sin(x), (x,0,10), plot_points=20, linestyle='', marker='.')

    The marker can be a TeX symbol as well::

        sage: plot(sin(x), (x,0,10), plot_points=20, linestyle='', marker=r'$\checkmark$')

    Sage currently ignores points that cannot be evaluated

    ::

        sage: set_verbose(-1)
        sage: plot(-x*log(x), (x,0,1))  # this works fine since the failed endpoint is just skipped.
        sage: set_verbose(0)

    This prints out a warning and plots where it can (we turn off the
    warning by setting the verbose mode temporarily to -1.)

    ::

        sage: set_verbose(-1)
        sage: plot(x^(1/3), (x,-1,1))
        sage: set_verbose(0)

    To plot the negative real cube root, use something like the following::

        sage: plot(lambda x : RR(x).nth_root(3), (x,-1, 1))

    We can detect the poles of a function::

        sage: plot(gamma, (-3, 4), detect_poles = True).show(ymin = -5, ymax = 5)

    We draw the Gamma-Function with its poles highlighted::

        sage: plot(gamma, (-3, 4), detect_poles = 'show').show(ymin = -5, ymax = 5)

    The basic options for filling a plot::

        sage: p1 = plot(sin(x), -pi, pi, fill = 'axis')
        sage: p2 = plot(sin(x), -pi, pi, fill = 'min')
        sage: p3 = plot(sin(x), -pi, pi, fill = 'max')
        sage: p4 = plot(sin(x), -pi, pi, fill = 0.5)
        sage: graphics_array([[p1, p2], [p3, p4]]).show(frame=True, axes=False)

        sage: plot([sin(x), cos(2*x)*sin(4*x)], -pi, pi, fill = {0: 1}, fillcolor = 'red', fillalpha = 1)

    A example about the growth of prime numbers::

        sage: plot(1.13*log(x), 1, 100, fill = lambda x: nth_prime(x)/floor(x), fillcolor = 'red')

    Fill the area between a function and its asymptote::

        sage: f = (2*x^3+2*x-1)/((x-2)*(x+1))
        sage: plot([f, 2*x+2], -7,7, fill = {0: [1]}, fillcolor='#ccc').show(ymin=-20, ymax=20)

    Fill the area between a list of functions and the x-axis::

        sage: def b(n): return lambda x: bessel_J(n, x)
        sage: plot([b(n) for n in [1..5]], 0, 20, fill = 'axis')

    Note that to fill between the ith and jth functions, you
    must use dictionary key-value pairs i:[j]; key-value pairs
    like i:j will fill between the ith function and the line y=j::

        sage: def b(n): return lambda x: bessel_J(n, x) + 0.5*(n-1)
        sage: plot([b(c) for c in [1..5]], 0, 40, fill = dict([(i, [i+1]) for i in [0..3]]))
        sage: plot([b(c) for c in [1..5]], 0, 40, fill = dict([(i, i+1) for i in [0..3]]))

    Extra options will get passed on to show(), as long as they are valid::

        sage: plot(sin(x^2), (x, -3, 3), axes_labels=['$x$','$y$']) # These labels will be nicely typeset
        sage: plot(sin(x^2), (x, -3, 3), axes_labels=['x','y']) # These will not

    ::

        sage: plot(sin(x^2), (x, -3, 3), figsize=[8,2])
        sage: plot(sin(x^2), (x, -3, 3)).show(figsize=[8,2]) # These are equivalent

    This includes options for custom ticks and formatting.  See documentation
    for :meth:`show` for more details.

    ::

        sage: plot(sin(pi*x), (x, -8, 8), ticks=[[-7,-3,0,3,7],[-1/2,0,1/2]])
        sage: plot(2*x+1,(x,0,5),ticks=[[0,1,e,pi,sqrt(20)],2],tick_formatter="latex")

    This is particularly useful when setting custom ticks in multiples of `pi`.

    ::

        sage: plot(sin(x),(x,0,2*pi),ticks=pi/3,tick_formatter=pi)

    A example with excluded values::

        sage: plot(floor(x), (x, 1, 10), exclude = [1..10])

    We exclude all points where prime_pi makes a jump::

        sage: jumps = [n for n in [1..100] if prime_pi(n) != prime_pi(n-1)]
        sage: plot(lambda x: prime_pi(x), (x, 1, 100), exclude = jumps)

    Excluded points can also be given by an equation::

        sage: g(x) = x^2-2*x-2
        sage: plot(1/g(x), (x, -3, 4), exclude = g(x) == 0, ymin = -5, ymax = 5)

    ``exclude`` and ``detect_poles`` can be used together::
        sage: f(x) = (floor(x)+0.5) / (1-(x-0.5)^2)
        sage: plot(f, (x, -3.5, 3.5), detect_poles = 'show', exclude = [-3..3], ymin = -5, ymax = 5)

    TESTS:

    We do not randomize the endpoints::

        sage: p = plot(x, (x,-1,1))
        sage: p[0].xdata[0] == -1
        True
        sage: p[0].xdata[-1] == 1
        True

    We check to make sure that the x/y min/max data get set correctly
    when there are multiple functions.

    ::

        sage: d = plot([sin(x), cos(x)], 100, 120).get_minmax_data()
        sage: d['xmin']
        100.0
        sage: d['xmax']
        120.0

    We check various combinations of tuples and functions, ending with
    tests that lambda functions work properly with explicit variable
    declaration, without a tuple.

    ::

        sage: p = plot(lambda x: x,(x,-1,1))
        sage: p = plot(lambda x: x,-1,1)
        sage: p = plot(x,x,-1,1)
        sage: p = plot(x,-1,1)
        sage: p = plot(x^2,x,-1,1)
        sage: p = plot(x^2,xmin=-1,xmax=2)
        sage: p = plot(lambda x: x,x,-1,1)
        sage: p = plot(lambda x: x^2,x,-1,1)
        sage: p = plot(lambda x: 1/x,x,-1,1)
        sage: f(x) = sin(x+3)-.1*x^3
        sage: p = plot(lambda x: f(x),x,-1,1)

    We check to handle cases where the function gets evaluated at a
    point which causes an 'inf' or '-inf' result to be produced.

    ::

        sage: p = plot(1/x, 0, 1)
        sage: p = plot(-1/x, 0, 1)

    Bad options now give better errors::

        sage: P = plot(sin(1/x), (x,-1,3), foo=10)
        Traceback (most recent call last):
        ...
        RuntimeError: Error in line(): option 'foo' not valid.
        sage: P = plot(x, (x,1,1)) # trac ticket #11753
        Traceback (most recent call last):
        ...
        ValueError: plot start point and end point must be different

    We test that we can plot `f(x)=x` (see Trac 10246)::

        sage: f(x)=x; f
        x |--> x
        sage: plot(f,(x,-1,1))
    """
    G_kwds = Graphics._extract_kwds_for_show(kwds, ignore=['xmin', 'xmax'])

    original_opts = kwds.pop('__original_opts', {})
    do_show = kwds.pop('show',False)

    from sage.structure.element import is_Vector
    if kwds.get('parametric',False) and is_Vector(funcs):
        funcs = tuple(funcs)


    if hasattr(funcs, 'plot'):
        G = funcs.plot(*args, **original_opts)
    # if we are using the generic plotting method
    else:
        n = len(args)
        # if there are no extra args, try to get xmin,xmax from
        # keyword arguments or pick some silly default
        if n == 0:
            xmin = kwds.pop('xmin', -1)
            xmax = kwds.pop('xmax', 1)
            G = _plot(funcs, (xmin, xmax), **kwds)

        # if there is one extra arg, then it had better be a tuple
        elif n == 1:
            G = _plot(funcs, *args, **kwds)
        elif n == 2:
        # if there are two extra args, then pull them out and pass them as a tuple
            xmin = args[0]
            xmax = args[1]
            args = args[2:]
            G = _plot(funcs, (xmin, xmax), *args, **kwds)
        elif n == 3:
        # if there are three extra args, then pull them out and pass them as a tuple
            var = args[0]
            xmin = args[1]
            xmax = args[2]
            args = args[3:]
            G = _plot(funcs, (var, xmin, xmax), *args, **kwds)
        elif ('xmin' in kwds) or ('xmax' in kwds):
            xmin = kwds.pop('xmin', -1)
            xmax = kwds.pop('xmax', 1)
            G = _plot(funcs, (xmin, xmax), *args, **kwds)
            pass
        else:
            sage.misc.misc.verbose("there were %s extra arguments (besides %s)" % (n, funcs), level=0)

    G._set_extra_kwds(G_kwds)
    if do_show:
        G.show()
    return G


def _plot(funcs, xrange, parametric=False,
              polar=False, fill=False, label='', randomize=True, **options):

    from sage.plot.misc import setup_for_eval_on_grid
    if funcs == []:
        return Graphics()
    funcs, ranges = setup_for_eval_on_grid(funcs, [xrange], options['plot_points'])
    xmin, xmax, delta = ranges[0]
    xrange=ranges[0][:2]
    #parametric_plot will be a list or tuple of two functions (f,g)
    #and will plotted as (f(x), g(x)) for all x in the given range
    if parametric:
        f, g = funcs
    #or we have only a single function to be plotted:
    else:
        f = funcs

    #check to see if funcs is a list of functions that will
    #be all plotted together.
    if isinstance(funcs, (list, tuple)) and not parametric:
        rainbow_colors = rainbow(len(funcs))

        G = Graphics()
        for i, h in enumerate(funcs):
            if isinstance(fill, dict):
                if i in fill:
                    fill_entry = fill[i]
                    if isinstance(fill_entry, list):
                        if fill_entry[0] < len(funcs):
                            fill_temp = funcs[fill_entry[0]]
                        else:
                            fill_temp = None
                    else:
                        fill_temp = fill_entry
                else:
                    fill_temp = None
            else:
                fill_temp = fill

            options_temp = options.copy()
            fillcolor_temp = options_temp.pop('fillcolor', 'automatic')
            legend_label=options_temp.pop('legend_label', None) # legend_label popped so the label isn't repeated for nothing
            if fillcolor_temp == 'automatic':
                fillcolor_temp = rainbow_colors[i]

            G += plot(h, xrange, polar = polar, fill = fill_temp, \
                      fillcolor = fillcolor_temp, **options_temp)
        return G

    adaptive_tolerance = options.pop('adaptive_tolerance')
    adaptive_recursion = options.pop('adaptive_recursion')
    plot_points = int(options.pop('plot_points'))

    exclude = options.pop('exclude')
    if exclude is not None:
        from sage.symbolic.expression import Expression
        if isinstance(exclude, Expression) and exclude.is_relational() == True:
            if len(exclude.variables()) > 1:
                raise ValueError('exclude has to be an equation of only one variable')
            v = exclude.variables()[0]
            points = [e.right() for e in exclude.solve(v) if e.left() == v and (v not in e.right().variables())]
            # We are only interested in real solutions
            exclude = []
            for x in points:
                try:
                    exclude.append(float(x))
                except TypeError:
                    pass

        if isinstance(exclude, (list, tuple)):
            exclude = sorted(exclude)
            # We make sure that points plot points close to the excluded points are computed
            epsilon = 0.001*(xmax - xmin)
            initial_points = reduce(lambda a,b: a+b, [[x - epsilon, x + epsilon] for x in exclude], [])
            data = generate_plot_points(f, xrange, plot_points, adaptive_tolerance, adaptive_recursion, randomize, initial_points)
        else:
            raise ValueError('exclude needs to be a list of numbers or an equation')

        if exclude == []:
            exclude = None
    else:
        data = generate_plot_points(f, xrange, plot_points, adaptive_tolerance, adaptive_recursion, randomize)

    if parametric:
        # We need the original x-values to be able to exclude points in parametric plots
        exclude_data = data
        data = [(fdata, g(x)) for x, fdata in data]

    G = Graphics()

    fillcolor = options.pop('fillcolor', 'automatic')
    fillalpha = options.pop('fillalpha', 0.5)

    # TODO: Use matplotlib's fill and fill_between commands.
    if fill is not False and fill is not None:
        if parametric:
            filldata = data
        else:
            if fill == 'axis' or fill is True:
                base_level = 0
            elif fill == 'min':
                base_level = min(t[1] for t in data)
            elif fill == 'max':
                base_level = max(t[1] for t in data)
            elif hasattr(fill, '__call__'):
                if fill == max or fill == min:
                    if fill == max:
                        fstr = 'max'
                    else:
                        fstr = 'min'
                    msg = "WARNING: You use the built-in function %s for filling. You probably wanted the string '%s'." % (fstr, fstr)
                    sage.misc.misc.verbose(msg, level=0)
                if not is_fast_float(fill):
                    fill_f = fast_float(fill, expect_one_var=True)
                else:
                    fill_f = fill

                filldata = generate_plot_points(fill_f, xrange, plot_points, adaptive_tolerance, \
                                                adaptive_recursion, randomize)
                filldata.reverse()
                filldata += data
            else:
                try:
                    base_level = float(fill)
                except TypeError:
                    base_level = 0

            if not hasattr(fill, '__call__') and polar:
                filldata = generate_plot_points(lambda x: base_level, xrange, plot_points, adaptive_tolerance, \
                                                adaptive_recursion, randomize)
                filldata.reverse()
                filldata += data
            if not hasattr(fill, '__call__') and not polar:
                filldata = [(data[0][0], base_level)] + data + [(data[-1][0], base_level)]

        if fillcolor == 'automatic':
            fillcolor = (0.5, 0.5, 0.5)
        fill_options = {}
        fill_options['rgbcolor'] = fillcolor
        fill_options['alpha'] = fillalpha
        fill_options['thickness'] = 0
        if polar:
            filldata = [(y*cos(x), y*sin(x)) for x, y in filldata]
        G += polygon(filldata, **fill_options)

    # We need the original data to be able to exclude points in polar plots
    if not parametric:
        exclude_data = data
    if polar:
        data = [(y*cos(x), y*sin(x)) for x, y in data]

    from sage.plot.all import line, text

    detect_poles = options.pop('detect_poles', False)
    if exclude is not None or detect_poles != False:
        start_index = 0
        # setup for pole detection
        from sage.rings.all import RDF
        epsilon = 0.0001
        pole_options = {}
        pole_options['linestyle'] = '--'
        pole_options['thickness'] = 1
        pole_options['rgbcolor'] = '#ccc'

        # setup for exclusion points
        exclusion_point = 0
        if exclude is not None:
            exclude.reverse()
            exclusion_point = exclude.pop()

        for i in range(len(data)-1):
            x0, y0 = exclude_data[i]
            x1, y1 = exclude_data[i+1]
            # detect poles
            if (not (polar or parametric)) and detect_poles != False \
               and ((y1 > 0 and y0 < 0) or (y1 < 0 and y0 > 0)):
                # calculate the slope of the line segment
                dy = abs(y1-y0)
                dx = x1 - x0
                alpha = (RDF(dy)/RDF(dx)).arctan()
                if alpha >= RDF(pi/2) - epsilon:
                    G += line(data[start_index:i], **options)
                    if detect_poles == 'show':
                        # draw a vertical asymptote
                        G += line([(x0, y0), (x1, y1)], **pole_options)
                    start_index = i+2

            # exclude points
            if exclude is not None and (x0 <= exclusion_point <= x1):
                G += line(data[start_index:i], **options)
                start_index = i + 2
                try:
                    exclusion_point = exclude.pop()
                except IndexError:
                    # all excluded points were considered
                    exclude = None

        G += line(data[start_index:], **options)
    else:
        G += line(data, **options)

    # Label?
    if label:
        sage.misc.misc.deprecation("Consider using legend_label instead")
        label = '  '+str(label)
        G += text(label, data[-1], horizontal_alignment='left',
                  vertical_alignment='center')

    return G




########## misc functions ###################

@options(aspect_ratio=1.0)
def parametric_plot(funcs, *args, **kwargs):
    r"""
    Plot a parametric curve or surface in 2d or 3d.

    :func:`parametric_plot` takes two or three functions as a
    list or a tuple and makes a plot with the first function giving the
    `x` coordinates, the second function giving the `y`
    coordinates, and the third function (if present) giving the
    `z` coordinates.

    In the 2d case, :func:`parametric_plot` is equivalent to the :func:`plot` command
    with the option ``parametric=True``.  In the 3d case, :func:`parametric_plot`
    is equivalent to :func:`~sage.plot.plot3d.parametric_plot3d.parametric_plot3d`.
    See each of these functions for more help and examples.

    INPUT:


    -  ``funcs`` - 2 or 3-tuple of functions, or a vector of dimension 2 or 3.

    -  ``other options`` - passed to :func:`plot` or :func:`~sage.plot.plot3d.parametric_plot3d.parametric_plot3d`


    EXAMPLES: We draw some 2d parametric plots.  Note that the default aspect ratio
    is 1, so that circles look like circles. ::

        sage: t = var('t')
        sage: parametric_plot( (cos(t), sin(t)), (t, 0, 2*pi))

    ::

        sage: parametric_plot( (sin(t), sin(2*t)), (t, 0, 2*pi), color=hue(0.6) )

    ::

        sage: parametric_plot((1, t), (t, 0, 4))

    Note that in parametric_plot, there is only fill or no fill.

    ::

        sage: parametric_plot((t, t^2), (t, -4, 4), fill = True)

    A filled Hypotrochoid::

        sage: parametric_plot([cos(x) + 2 * cos(x/4), sin(x) - 2 * sin(x/4)], (x,0, 8*pi), fill = True)

        sage: parametric_plot( (5*cos(x), 5*sin(x), x), (x,-12, 12), plot_points=150, color="red")

        sage: y=var('y')
        sage: parametric_plot( (5*cos(x), x*y, cos(x*y)), (x, -4,4), (y,-4,4))

        sage: t=var('t')
        sage: parametric_plot( vector((sin(t), sin(2*t))), (t, 0, 2*pi), color='green')
        sage: parametric_plot( vector([t, t+1, t^2]), (t, 0, 1))

    TESTS::

        sage: parametric_plot((x, t^2), (x, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((1, x+t), (x, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((-t, x+t), (x, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((1, x+t, y), (x, -4, 4), (t, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((1, x, y), 0, 4)
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges
    """
    num_ranges=0
    for i in args:
        if isinstance(i, (list, tuple)):
            num_ranges+=1
        else:
            break

    if num_ranges==0 and len(args)>=2:
        from sage.misc.misc import deprecation
        deprecation("variable ranges to parametric_plot must be given as tuples, like (2,4) or (t,2,3)")
        args=tuple(args)
        num_ranges=1

    num_funcs = len(funcs)

    num_vars=len(sage.plot.misc.unify_arguments(funcs)[0])
    if num_vars>num_ranges:
        raise ValueError, "there are more variables than variable ranges"

    if num_funcs == 2 and num_ranges == 1:
        kwargs['parametric'] = True
        return plot(funcs, *args, **kwargs)
    elif (num_funcs == 3 and num_ranges <= 2):
        return sage.plot.plot3d.parametric_plot3d.parametric_plot3d(funcs, *args, **kwargs)
    else:
        raise ValueError, "the number of functions and the number of variable ranges is not a supported combination for a 2d or 3d parametric plots"

@options(aspect_ratio=1.0)
def polar_plot(funcs, *args, **kwds):
    r"""
    ``polar_plot`` takes a single function or a list or
    tuple of functions and plots them with polar coordinates in the given
    domain.

    This function is equivalent to the :func:`plot` command with the options
    ``polar=True`` and ``aspect_ratio=1``. For more help on options,
    see the documentation for :func:`plot`.

    INPUT:

    - ``funcs`` - a function
    - other options are passed to plot

    EXAMPLES:

    Here is a blue 8-leaved petal::

        sage: polar_plot(sin(5*x)^2, (x, 0, 2*pi), color='blue')

    A red figure-8::

        sage: polar_plot(abs(sqrt(1 - sin(x)^2)), (x, 0, 2*pi), color='red')

    A green limacon of Pascal::

        sage: polar_plot(2 + 2*cos(x), (x, 0, 2*pi), color=hue(0.3))

    Several polar plots::

        sage: polar_plot([2*sin(x), 2*cos(x)], (x, 0, 2*pi))

    A filled spiral::

        sage: polar_plot(sqrt, 0, 2 * pi, fill = True)

    Fill the area between two functions::

        sage: polar_plot(cos(4*x) + 1.5, 0, 2*pi, fill=0.5 * cos(4*x) + 2.5, fillcolor='orange')

    Fill the area between several spirals::

        sage: polar_plot([(1.2+k*0.2)*log(x) for k in range(6)], 1, 3 * pi, fill = {0: [1], 2: [3], 4: [5]})

    Exclude points at discontinuities::

        sage: polar_plot(log(floor(x)), (x, 1, 4*pi), exclude = [1..12])

    """
    kwds['polar']=True
    return plot(funcs, *args, **kwds)

@options(aspect_ratio='automatic')
def list_plot(data, plotjoined=False, **kwargs):
    r"""
    ``list_plot`` takes either a list of numbers, a list of tuples,
    or a dictionary and plots the corresponding points.

    If given a list of numbers (that is, not a list of tuples or lists),
    ``list_plot`` forms a list of tuples `(i, x_i)` where `i` goes from
    0 to ``len(data)-1`` and `x_i` is the `i`-th data value, and puts
    points at those tuple values.

    ``list_plot`` will plot a list of complex numbers in the obvious
    way; any numbers for which
    :func:`CC()<sage.rings.complex_field.ComplexField>` makes sense will
    work.

    ``list_plot`` also takes a list of tuples `(x_i, y_i)` where `x_i`
    and `y_i` are the `i`-th values representing the `x`- and
    `y`-values, respectively.

    If given a dictionary, ``list_plot`` interprets the keys as
    `x`-values and the values as `y`-values.

    The ``plotjoined=True`` option tells ``list_plot`` to plot a line
    joining all the data.

    It is possible to pass empty dictionaries, lists, or tuples to
    list_plot. Doing so will plot nothing (returning an empty plot).

    EXAMPLES::

        sage: list_plot([i^2 for i in range(5)])

    Here are a bunch of random red points::

        sage: r = [(random(),random()) for _ in range(20)]
        sage: list_plot(r,color='red')

    This gives all the random points joined in a purple line::

        sage: list_plot(r, plotjoined=True, color='purple')

    Plot a list of complex numbers::

        sage: list_plot([1, I, pi + I/2, CC(.25, .25)])

        sage: list_plot([exp(I*theta) for theta in [0, .2..pi]])

    Note that if your list of complex numbers are all actually real,
    they get plotted as real values, so this

    ::

        sage: list_plot([CDF(1), CDF(1/2), CDF(1/3)])

    is the same as ``list_plot([1, 1/2, 1/3])`` -- it produces a plot of
    the points `(0,1)`, `(1,1/2)`, and `(2,1/3)`.

    If you have separate lists of `x` values and `y` values which you
    want to plot against each other, use the ``zip`` command to make a
    single list whose entries are pairs of `(x,y)` values, and feed
    the result into ``list_plot``::

        sage: x_coords = [cos(t)^3 for t in srange(0, 2*pi, 0.02)]
        sage: y_coords = [sin(t)^3 for t in srange(0, 2*pi, 0.02)]
        sage: list_plot(zip(x_coords, y_coords))

    If instead you try to pass the two lists as separate arguments,
    you will get an error message::

        sage: list_plot(x_coords, y_coords)
        Traceback (most recent call last):
        ...
        TypeError: The second argument 'plotjoined' should be boolean (True or False).  If you meant to plot two lists 'x' and 'y' against each other, use 'list_plot(zip(x,y))'.

    Dictionaries with numeric keys and values can be plotted::

        sage: list_plot({22: 3365, 27: 3295, 37: 3135, 42: 3020, 47: 2880, 52: 2735, 57: 2550})

    TESTS:

    We check to see that the x/y min/max data are set correctly.

    ::

        sage: d = list_plot([(100,100), (120, 120)]).get_minmax_data()
        sage: d['xmin']
        100.0
        sage: d['ymin']
        100.0
    """
    from sage.plot.all import line, point
    if data == {} or data == () or data == []:
        return Graphics()
    if isinstance(data, dict):
        if plotjoined:
            list_data = sorted(list(data.iteritems()))
        else:
            list_data = list(data.iteritems())
        return list_plot(list_data, plotjoined=plotjoined, **kwargs)
    if not isinstance(data[0], (list, tuple)):
        data = zip(range(len(data)), data)
    if isinstance(plotjoined, (list, tuple)):
        raise TypeError, "The second argument 'plotjoined' should be boolean (True or False).  If you meant to plot two lists 'x' and 'y' against each other, use 'list_plot(zip(x,y))'."
    try:
        if plotjoined:
            return line(data, **kwargs)
        else:
            return point(data, **kwargs)
    except (TypeError, IndexError):
        # Assume we have complex-valued input and plot real and imaginary parts.
        # Need to catch IndexError because if data is, say, [(0, 1), (1, I)],
        # point3d() throws an IndexError on the (0,1) before it ever
        # gets to (1, I).
        from sage.rings.complex_field import ComplexField
        CC = ComplexField()
        # if we get here, we already did "zip(range(len(data)), data)",
        # so look at z[1] in inner list
        data = [(z.real(), z.imag()) for z in [CC(z[1]) for z in data]]
        if plotjoined:
            return line(data, **kwargs)
        else:
            return point(data, **kwargs)

def to_float_list(v):
    """
    Given a list or tuple or iterable v, coerce each element of v to a
    float and make a list out of the result.

    EXAMPLES::

        sage: from sage.plot.plot import to_float_list
        sage: to_float_list([1,1/2,3])
        [1.0, 0.5, 3.0]
    """
    return [float(x) for x in v]

class GraphicsArray(SageObject):
    """
    GraphicsArray takes a (`m` x `n`) list of lists of
    graphics objects and plots them all on one canvas.
    """
    def __init__(self, array):
        if not isinstance(array, (list, tuple)):
            raise TypeError,"array (=%s) must be a list of lists of Graphics objects"%(array)
        array = list(array)
        self._glist = []
        self._rows = len(array)
        if self._rows > 0:
            if not isinstance(array[0], (list, tuple)):
                array = [array]
                self._rows = 1
            self._cols = len(array[0])
        else:
            self._cols = 0
        self._dims = self._rows*self._cols
        for row in array: #basically flatten the list
            if not isinstance(row, (list, tuple)) or len(row) != self._cols:
                raise TypeError,"array (=%s) must be a list of lists of Graphics objects"%(array)
            for g in row:
                if not isinstance(g, Graphics):
                    raise TypeError, "every element of array must be a Graphics object"
                self._glist.append(g)
        self._figsize = None

    def _repr_(self):
        if SHOW_DEFAULT:
            self.show()
            return ''
        else:
            return self.__str__()

    def __str__(self):
        return "Graphics Array of size %s x %s"%(self._rows, self._cols)

    def nrows(self):
        return self._rows

    def ncols(self):
        return self._cols

    def __getitem__(self, i):
        i = int(i)
        return self._glist[i]

    def __setitem__(self, i, g):
        i = int(i)
        self._glist[i] = g

    def __set_figsize__(self, list):
        m = int(list[0])
        n = int(list[1])
        self._figsize = [m,n]

    def __len__(self):
        return len(self._glist)

    def append(self, g):
        self._glist.append(g)

    def _render(self, filename, dpi=None, figsize=None, axes=None, **args):
        r"""
        ``render`` loops over all graphics objects in the array
        and adds them to the subplot.
        """
        #glist is a list of Graphics objects:
        glist = self._glist
        rows = self._rows
        cols = self._cols
        dims = self._dims
        #make a blank matplotlib Figure:
        from matplotlib.figure import Figure
        figure = Figure(figsize)
        global do_verify
        do_verify = True
        for i,g in zip(range(1, dims+1), glist):
            subplot = figure.add_subplot(rows, cols, i)
            g.matplotlib(filename, figure=figure, sub=subplot,
                         verify=do_verify, axes = axes, **args)
        g.save(filename, dpi=dpi, figure=figure, sub=subplot,
               verify=do_verify, axes = axes, **args)

    def save(self, filename=None, dpi=DEFAULT_DPI, figsize=None,
             axes = None, **args):
        """
        save the ``graphics_array`` to (for now) a png called
        'filename'.
        """
        if (figsize is not None): self.__set_figsize__(figsize)
        self._render(filename, dpi=dpi, figsize=self._figsize, axes = axes, **args)

    def show(self, filename=None, dpi=DEFAULT_DPI, figsize=None,
             axes = None, **args):
        r"""
        Show this graphics array using the default viewer.

        OPTIONAL INPUT:


        -  ``filename`` - (default: None) string

        -  ``dpi`` - dots per inch

        -  ``figsize`` - width or [width, height]

        -  ``axes`` - (default: True)

        -  ``fontsize`` - positive integer

        -  ``frame`` - (default: False) draw a frame around the
           image


        EXAMPLES: This draws a graphics array with four trig plots and no
        axes in any of the plots.

        ::

            sage: G = graphics_array([[plot(sin), plot(cos)], [plot(tan), plot(sec)]])
            sage: G.show(axes=False)
        """
        if (figsize is not None): self.__set_figsize__(figsize)
        if DOCTEST_MODE:
            self.save(DOCTEST_MODE_FILE,
                      dpi=dpi, figsize=self._figsize, axes = axes, **args)
            return
        if EMBEDDED_MODE:
            self.save(filename, dpi=dpi, figsize=self._figsize, axes = axes, **args)
            return
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self._render(filename, dpi=dpi, figsize=self._figsize, axes = axes, **args)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(
                         sage.misc.viewer.browser(), filename))


def reshape(v, n, m):
    G = Graphics()
    G.axes(False)
    if len(v) == 0:
        return [[G]*m]*n

    if not isinstance(v[0], Graphics):
        # a list of lists -- flatten it
        v = sum([list(x) for x in v], [])

    # Now v should be a single list.
    # First, make it have the right length.
    for i in xrange(n*m - len(v)):
        v.append(G)

    # Next, create a list of lists out of it.
    L = []
    k = 0
    for i in range(n):
        w = []
        for j in range(m):
            w.append(v[k])
            k += 1
        L.append(w)

    return L

def graphics_array(array, n=None, m=None):
    r"""
    ``graphics_array`` take a list of lists (or tuples) of
    graphics objects and plots them all on one canvas (single plot).

    INPUT:


    -  ``array`` - a list of lists or tuples

    -  ``n, m`` - (optional) integers - if n and m are
       given then the input array is flattened and turned into an n x m
       array, with blank graphics objects padded at the end, if
       necessary.


    EXAMPLE: Make some plots of `\sin` functions::

        sage: f(x) = sin(x)
        sage: g(x) = sin(2*x)
        sage: h(x) = sin(4*x)
        sage: p1 = plot(f,(-2*pi,2*pi),color=hue(0.5))
        sage: p2 = plot(g,(-2*pi,2*pi),color=hue(0.9))
        sage: p3 = parametric_plot((f,g),(0,2*pi),color=hue(0.6))
        sage: p4 = parametric_plot((f,h),(0,2*pi),color=hue(1.0))

    Now make a graphics array out of the plots; then you can type
    either: ``ga.show()`` or ``ga.save()``.

    ::

        sage: graphics_array(((p1,p2),(p3,p4)))

    Here we give only one row::

        sage: p1 = plot(sin,(-4,4))
        sage: p2 = plot(cos,(-4,4))
        sage: g = graphics_array([p1, p2]); print g
        Graphics Array of size 1 x 2
        sage: g.show()
    """
    if not n is None:
        # Flatten then reshape input
        n = int(n)
        m = int(m)
        array = reshape(array, n, m)
    return GraphicsArray(array)

def var_and_list_of_values(v, plot_points):
    """
    INPUT:


    -  ``v`` - (v0, v1) or (var, v0, v1); if the former
       return the range of values between v0 and v1 taking plot_points
       steps; if var is given, also return var.

    -  ``plot_points`` - integer = 2 (the endpoints)


    OUTPUT:


    -  ``var`` - a variable or None

    -  ``list`` - a list of floats


    EXAMPLES::

        sage: from sage.plot.plot import var_and_list_of_values
        sage: var_and_list_of_values((var('theta'), 2, 5),  5)
        doctest:...: DeprecationWarning: var_and_list_of_values is deprecated.  Please use sage.plot.misc.setup_for_eval_on_grid; note that that function has slightly different calling and return conventions which make it more generally applicable
        (theta, [2.0, 2.75, 3.5, 4.25, 5.0])
        sage: var_and_list_of_values((2, 5),  5)
        (None, [2.0, 2.75, 3.5, 4.25, 5.0])
        sage: var_and_list_of_values((var('theta'), 2, 5),  2)
        (theta, [2.0, 5.0])
        sage: var_and_list_of_values((2, 5),  2)
        (None, [2.0, 5.0])
    """
    from sage.misc.misc import deprecation
    deprecation("var_and_list_of_values is deprecated.  Please use sage.plot.misc.setup_for_eval_on_grid; note that that function has slightly different calling and return conventions which make it more generally applicable")
    plot_points = int(plot_points)
    if plot_points < 2:
        raise ValueError, "plot_points must be greater than 1"
    if not isinstance(v, (tuple, list)):
        raise TypeError, "v must be a tuple or list"
    if len(v) == 3:
        var = v[0]
        a, b = v[1], v[2]
    elif len(v) == 2:
        var = None
        a, b = v
    else:
        raise ValueError, "parametric value range must be a list or tuple of length 2 or 3."

    a = float(a)
    b = float(b)
    if plot_points == 2:
        return var, [a, b]
    else:
        step = (b-a)/float(plot_points-1)
        values = [a + step*i for i in xrange(plot_points)]
        return var, values



def setup_for_eval_on_grid(v, xrange, yrange, plot_points):
    """
    This function is deprecated.  Please use
    sage.plot.misc.setup_for_eval_on_grid instead.  Please note that
    that function has slightly different calling and return
    conventions which make it more generally applicable.

    INPUT:


    -  ``v`` - a list of functions

    -  ``xrange`` - 2 or 3 tuple (if 3, first is a
       variable)

    -  ``yrange`` - 2 or 3 tuple

    -  ``plot_points`` - a positive integer


    OUTPUT:


    -  ``g`` - tuple of fast callable functions

    -  ``xstep`` - step size in xdirection

    -  ``ystep`` - step size in ydirection

    -  ``xrange`` - tuple of 2 floats

    -  ``yrange`` - tuple of 2 floats


    EXAMPLES::

        sage: x,y = var('x,y')
        sage: sage.plot.plot.setup_for_eval_on_grid([x^2 + y^2], (x,0,5), (y,0,pi), 11)
        doctest:...: DeprecationWarning: sage.plot.plot.setup_for_eval_on_grid is deprecated.  Please use sage.plot.misc.setup_for_eval_on_grid; note that that function has slightly different calling and return conventions which make it more generally applicable
        ([<sage.ext... object at ...>],
         0.5,
         0.3141592653589793,
         (0.0, 5.0),
         (0.0, 3.141592653589793))

    We always plot at least two points; one at the beginning and one at the end of the ranges.

    ::

        sage: sage.plot.plot.setup_for_eval_on_grid([x^2+y^2], (x,0,1), (y,-1,1), 1)
        ([<sage.ext... object at ...>],
        1.0,
        2.0,
        (0.0, 1.0),
        (-1.0, 1.0))


    """
    from sage.misc.misc import deprecation
    deprecation("sage.plot.plot.setup_for_eval_on_grid is deprecated.  Please use sage.plot.misc.setup_for_eval_on_grid; note that that function has slightly different calling and return conventions which make it more generally applicable")

    from sage.plot.misc import setup_for_eval_on_grid as setup
    g, ranges=setup(v, [xrange, yrange], plot_points)
    return list(g), ranges[0][2], ranges[1][2], ranges[0][:2], ranges[1][:2]


def minmax_data(xdata, ydata, dict=False):
    """
    Returns the minimums and maximums of xdata and ydata.

    If dict is False, then minmax_data returns the tuple (xmin, xmax,
    ymin, ymax); otherwise, it returns a dictionary whose keys are
    'xmin', 'xmax', 'ymin', and 'ymax' and whose values are the
    corresponding values.

    EXAMPLES::

        sage: from sage.plot.plot import minmax_data
        sage: minmax_data([], [])
        (-1, 1, -1, 1)
        sage: minmax_data([-1, 2], [4, -3])
        (-1, 2, -3, 4)
        sage: d = minmax_data([-1, 2], [4, -3], dict=True)
        sage: list(sorted(d.items()))
        [('xmax', 2), ('xmin', -1), ('ymax', 4), ('ymin', -3)]
    """
    xmin = min(xdata) if len(xdata) > 0 else -1
    xmax = max(xdata) if len(xdata) > 0 else 1
    ymin = min(ydata) if len(ydata) > 0 else -1
    ymax = max(ydata) if len(ydata) > 0 else 1
    if dict:
        return {'xmin':xmin, 'xmax':xmax,
                'ymin':ymin, 'ymax':ymax}
    else:
        return xmin, xmax, ymin, ymax

def adaptive_refinement(f, p1, p2, adaptive_tolerance=0.01, adaptive_recursion=5, level=0):
    r"""
    The adaptive refinement algorithm for plotting a function f. See
    the docstring for plot for a description of the algorithm.

    INPUT:


    -  ``f`` - a function of one variable

    -  ``p1, p2`` - two points to refine between

    -  ``adaptive_recursion`` - (default: 5) how many
       levels of recursion to go before giving up when doing adaptive
       refinement. Setting this to 0 disables adaptive refinement.

    -  ``adaptive_tolerance`` - (default: 0.01) how large
       a relative difference should be before the adaptive refinement
       code considers it significant; see documentation for generate_plot_points
       for more information.  See the documentation for plot() for more
       information on how the adaptive refinement algorithm works.

    OUTPUT:


    -  ``list`` - a list of points to insert between p1 and
       p2 to get a better linear approximation between them


    TESTS::

        sage: from sage.plot.plot import adaptive_refinement
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01, adaptive_recursion=0)
        []
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01)
        [(0.125*pi, 0.3826834323650898), (0.1875*pi, 0.5555702330196022), (0.25*pi, 0.7071067811865475), (0.3125*pi, 0.8314696123025452), (0.375*pi, 0.9238795325112867), (0.4375*pi, 0.9807852804032304), (0.5*pi, 1.0), (0.5625*pi, 0.9807852804032304), (0.625*pi, 0.9238795325112867), (0.6875*pi, 0.8314696123025455), (0.75*pi, 0.7071067811865476), (0.8125*pi, 0.5555702330196022), (0.875*pi, 0.3826834323650899)]

    This shows that lowering adaptive_tolerance and raising
    adaptive_recursion both increase the number of subdivision
    points, though which one creates more points is heavily
    dependent upon the function being plotted.

    ::

        sage: x = var('x')
        sage: f(x) = sin(1/x)
        sage: n1 = len(adaptive_refinement(f, (0,0), (pi,0), adaptive_tolerance=0.01)); n1
        15
        sage: n2 = len(adaptive_refinement(f, (0,0), (pi,0), adaptive_recursion=10, adaptive_tolerance=0.01)); n2
        79
        sage: n3 = len(adaptive_refinement(f, (0,0), (pi,0), adaptive_tolerance=0.001)); n3
        26
    """
    if level >= adaptive_recursion:
        return []

    x = (p1[0] + p2[0])/2.0
    msg = ''

    try:
        y = float(f(x))
        if str(y) in ['nan', 'NaN', 'inf', '-inf']:
            sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, x),1)
            # give up for this branch
            return []

    except (ZeroDivisionError, TypeError, ValueError, OverflowError), msg:
        sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, x), 1)
        # give up for this branch
        return []

    # this distance calculation is not perfect.
    if abs((p1[1] + p2[1])/2.0 - y) > adaptive_tolerance:
        return adaptive_refinement(f, p1, (x, y),
                    adaptive_tolerance=adaptive_tolerance,
                    adaptive_recursion=adaptive_recursion,
                    level=level+1) \
                    + [(x, y)] + \
            adaptive_refinement(f, (x, y), p2,
                    adaptive_tolerance=adaptive_tolerance,
                    adaptive_recursion=adaptive_recursion,
                    level=level+1)
    else:
        return []

def generate_plot_points(f, xrange, plot_points=5, adaptive_tolerance=0.01, adaptive_recursion=5, randomize = True, initial_points = None):
    r"""
    Calculate plot points for a function f in the interval xrange.  The
    adaptive refinement algorithm is also automatically invoked with a
    *relative* adaptive tolerance of adaptive_tolerance; see below.

    INPUT:

    - ``f`` - a function of one variable

    - ``p1, p2`` - two points to refine between

    - ``plot_points`` - (default: 5) the minimal number of plot points. (Note
      however that in any actual plot a number is passed to this, with default
      value 200.)

    - ``adaptive_recursion`` - (default: 5) how many levels of recursion to go
      before giving up when doing adaptive refinement.  Setting this to 0
      disables adaptive refinement.

    - ``adaptive_tolerance`` - (default: 0.01) how large the relative difference
      should be before the adaptive refinement code considers it significant.  If
      the actual difference is greater than adaptive_tolerance*delta, where delta
      is the initial subinterval size for the given xrange and plot_points, then
      the algorithm will consider it significant.

    - ``initial_points`` - (default: None) a list of points that should be evaluated.

    OUTPUT:

    - a list of points (x, f(x)) in the interval xrange, which approximate
      the function f.

    TESTS::

        sage: from sage.plot.plot import generate_plot_points
        sage: generate_plot_points(sin, (0, pi), plot_points=2, adaptive_recursion=0)
        [(0.0, 0.0), (3.141592653589793, 1.2246...e-16)]

        sage: from sage.plot.plot import generate_plot_points
        sage: generate_plot_points(lambda x: x^2, (0, 6), plot_points=2, adaptive_recursion=0, initial_points = [1,2,3])
        [(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0), (6.0, 36.0)]

        sage: generate_plot_points(sin(x).function(x), (-pi, pi), randomize=False)
        [(-3.141592653589793, -1.2246...e-16), (-2.748893571891069,
        -0.3826834323650899), (-2.356194490192345, -0.707106781186547...),
        (-2.1598449493429825, -0.831469612302545...), (-1.9634954084936207,
        -0.9238795325112867), (-1.7671458676442586, -0.9807852804032304),
        (-1.5707963267948966, -1.0), (-1.3744467859455345,
        -0.9807852804032304), (-1.1780972450961724, -0.9238795325112867),
        (-0.9817477042468103, -0.831469612302545...), (-0.7853981633974483,
        -0.707106781186547...), (-0.39269908169872414, -0.3826834323650898),
        (0.0, 0.0), (0.39269908169872414, 0.3826834323650898),
        (0.7853981633974483, 0.707106781186547...), (0.9817477042468103,
        0.831469612302545...), (1.1780972450961724, 0.9238795325112867),
        (1.3744467859455345, 0.9807852804032304), (1.5707963267948966, 1.0),
        (1.7671458676442586, 0.9807852804032304), (1.9634954084936207,
        0.9238795325112867), (2.1598449493429825, 0.831469612302545...),
        (2.356194490192345, 0.707106781186547...), (2.748893571891069,
        0.3826834323650899), (3.141592653589793, 1.2246...e-16)]

    This shows that lowering adaptive_tolerance and raising
    adaptive_recursion both increase the number of subdivision points.
    (Note that which creates more points is heavily dependent on the
    particular function plotted.)

    ::

        sage: x = var('x')
        sage: f(x) = sin(1/x)
        sage: [len(generate_plot_points(f, (-pi, pi), plot_points=16, adaptive_tolerance=i, randomize=False)) for i in [0.01, 0.001, 0.0001]]
        [97, 161, 275]

        sage: [len(generate_plot_points(f, (-pi, pi), plot_points=16, adaptive_recursion=i, randomize=False)) for i in [5, 10, 15]]
        [97, 499, 2681]
    """
    from sage.plot.misc import setup_for_eval_on_grid
    ignore, ranges = setup_for_eval_on_grid([], [xrange], plot_points)
    xmin, xmax, delta = ranges[0]
    data = srange(*ranges[0], include_endpoint=True)

    random = current_randstate().python_random().random

    for i in range(len(data)):
        xi = data[i]
        # Slightly randomize the interior sample points if
        # randomize is true
        if randomize and i > 0 and i < plot_points-1:
            xi += delta*(random() - 0.5)
            data[i] = xi

    # add initial points
    if isinstance(initial_points, list):
        data = sorted(data + initial_points)

    exceptions = 0; msg=''
    exception_indices = []
    for i in range(len(data)):
        xi = data[i]

        try:
            data[i] = (float(xi), float(f(xi)))
            if str(data[i][1]) in ['nan', 'NaN', 'inf', '-inf']:
                sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, xi),1)
                exceptions += 1
                exception_indices.append(i)

        except (ArithmeticError, TypeError, ValueError), msg:
            sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, xi),1)

            if i == 0: # Given an error for left endpoint, try to move it in slightly
                for j in range(1, 99):
                    xj = xi + delta*j/100.0
                    try:
                        data[i] = (float(xj), float(f(xj)))
                        # nan != nan
                        if data[i][1] != data[i][1]:
                            continue
                        break
                    except (ArithmeticError, TypeError, ValueError), msg:
                        pass
                else:
                    exceptions += 1
                    exception_indices.append(i)

            elif i == plot_points-1: # Given an error for right endpoint, try to move it in slightly
                for j in range(1, 99):
                    xj = xi - delta*j/100.0
                    try:
                        data[i] = (float(xj), float(f(xj)))
                        # nan != nan
                        if data[i][1] != data[i][1]:
                            continue
                        break
                    except (ArithmeticError, TypeError, ValueError), msg:
                        pass
                else:
                    exceptions += 1
                    exception_indices.append(i)
            else:
                exceptions += 1
                exception_indices.append(i)

    data = [data[i] for i in range(len(data)) if i not in exception_indices]

    # calls adaptive refinement
    i, j = 0, 0
    adaptive_tolerance = delta * float(adaptive_tolerance)
    adaptive_recursion = int(adaptive_recursion)

    while i < len(data) - 1:
       for p in adaptive_refinement(f, data[i], data[i+1],
                                     adaptive_tolerance=adaptive_tolerance,
                                     adaptive_recursion=adaptive_recursion):
            data.insert(i+1, p)
            i += 1
       i += 1

    if (len(data) == 0 and exceptions > 0) or exceptions > 10:
        sage.misc.misc.verbose("WARNING: When plotting, failed to evaluate function at %s points."%exceptions, level=0)
        sage.misc.misc.verbose("Last error message: '%s'"%msg, level=0)

    return data

#Lovely cruft to keep the pickle jar working
from line import line, line2d, Line as GraphicPrimitive_Line
from arrow import arrow, Arrow as GraphicPrimitive_Arrow
from bar_chart import bar_chart, BarChart as GraphicPrimitive_BarChart
from disk import disk, Disk as GraphicPrimitive_Disk
from point import point, points, point2d, Point as GraphicPrimitive_Point
from matrix_plot import matrix_plot, MatrixPlot as GraphicPrimitive_MatrixPlot
from plot_field import plot_vector_field, plot_slope_field, PlotField as GraphicPrimitive_PlotField
from text import text, Text as GraphicPrimitive_Text
from polygon import polygon, Polygon as GraphicPrimitive_Polygon
from circle import circle, Circle as GraphicPrimtive_Circle
from contour_plot import contour_plot, implicit_plot, ContourPlot as GraphicPrimitive_ContourPlot

