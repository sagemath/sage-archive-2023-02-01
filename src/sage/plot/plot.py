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

-  The following log plotting functions:

   - :func:`plot_loglog`

   - :func:`plot_semilogx` and :func:`plot_semilogy`

   - :func:`list_plot_loglog`

   - :func:`list_plot_semilogx` and :func:`list_plot_semilogy`


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

You can even have custom tick labels along with custom positioning. ::

    sage: plot(x**2, (x,0,3), ticks=[[1,2.5],pi/2], tick_formatter=[["$x_1$","$x_2$"],pi]) # long time

We construct a plot involving several graphics objects::

    sage: G = plot(cos(x), (x, -5, 5), thickness=5, color='green', title='A plot')
    sage: P = polygon([[1,2], [5,6], [5,0]], color='red')
    sage: G + P

Next we construct the reflection of the above polygon about the
`y`-axis by iterating over the list of first-coordinates of
the first graphic element of ``P`` (which is the actual
Polygon; note that ``P`` is a Graphics object, which consists
of a single polygon)::

    sage: Q = polygon([(-x,y) for x,y in P[0]], color='blue')
    sage: Q   # show it

We combine together different graphics objects using "+"::

    sage: H = G + P + Q
    sage: print H
    Graphics object consisting of 3 graphics primitives
    sage: type(H)
    <class 'sage.plot.graphics.Graphics'>
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

We can add a title to a graph::

    sage: x = var('x')
    sage: plot(x^2, (x,-2,2), title='A plot of $x^2$')

We can set the position of the title::

    sage: plot(x^2, (-2,2), title='Plot of $x^2$', title_pos=(0.5,-0.05))

We plot the Riemann zeta function along the critical line and see
the first few zeros::

    sage: i = CDF.0      # define i this way for maximum speed.
    sage: p1 = plot(lambda t: arg(zeta(0.5+t*i)), 1,27,rgbcolor=(0.8,0,0))
    sage: p2 = plot(lambda t: abs(zeta(0.5+t*i)), 1,27,color=hue(0.7))
    sage: print p1 + p2
    Graphics object consisting of 2 graphics primitives
    sage: p1 + p2    # display it

Many concentric circles shrinking toward the origin::

    sage: show(sum(circle((i,0), i, hue=sin(i/10)) for i in [10,9.9,..,0])) # long time

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
    sage: (g1+g2).show(ticks=pi/6, tick_formatter=pi)  # long time # show their sum, nicely formatted

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
:trac:`2900` is fixed (in Matplotlib).

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

    sage: os.system("sage -c \"if 'matplotlib' in sys.modules: sys.exit(1)\"") # long time
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

- Jeroen Demeyer (2012-04-19): move parts of this file to graphics.py (:trac:`12857`)

"""
#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>
#       Copyright (C) 2006-2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Jason Grout
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os
from functools import reduce

## IMPORTANT: Do *not* import matplotlib at module scope.  It takes a
## surprisingly long time to initialize itself.  It's better if it is
## imported in functions, so it only gets started if it is actually
## going to be used.

#DEFAULT_FIGSIZE=(6, 3.70820393249937)
EMBEDDED_MODE = False
import sage.misc.misc
from sage.misc.misc import srange

from sage.misc.randstate import current_randstate #for plot adaptive refinement
from math import sin, cos, pi #for polar_plot

from sage.ext.fast_eval import fast_float, fast_float_constant, is_fast_float

from sage.misc.decorators import options, rename_keyword

from graphics import Graphics, GraphicsArray

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
        raise ValueError("points must have 2 coordinates in a 2d line")


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
    graphics primitives like :class:`~sage.plot.arrow.Arrow`;
    type ``sage.plot.plot?`` for a current list.

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

    - ``base`` - (default: 10) the base of the logarithm if
      a logarithmic scale is set. This must be greater than 1. The base
      can be also given as a list or tuple ``(basex, basey)``.
      ``basex`` sets the base of the logarithm along the horizontal
      axis and ``basey`` sets the base along the vertical axis.

    - ``scale`` -- (default: ``"linear"``) string. The scale of the axes.
      Possible values are ``"linear"``, ``"loglog"``, ``"semilogx"``,
      ``"semilogy"``.

      The scale can be also be given as single argument that is a list
      or tuple ``(scale, base)`` or ``(scale, basex, basey)``.

      The ``"loglog"`` scale sets both the horizontal and vertical axes to
      logarithmic scale. The ``"semilogx"`` scale sets the horizontal axis
      to logarithmic scale. The ``"semilogy"`` scale sets the vertical axis
      to logarithmic scale. The ``"linear"`` scale is the default value
      when :class:`~sage.plot.graphics.Graphics` is initialized.

    - ``xmin`` - starting x value

    - ``xmax`` - ending x value

    - ``ymin`` - starting y value in the rendered figure

    - ``ymax`` - ending y value in the rendered figure

    - ``color`` - an RGB tuple (r,g,b) with each of r,g,b between 0 and 1,
      or a color name as a string (e.g., 'purple'), or an HTML color
      such as '#aaff0b'.

    - ``detect_poles`` - (Default: False) If set to True poles are detected.
      If set to "show" vertical asymptotes are drawn.

    - ``legend_color`` - the color of the text for this item in the legend

    - ``legend_label`` - the label for this item in the legend

    .. note::

        - If the ``scale`` is ``"linear"``, then irrespective of what
          ``base`` is set to, it will default to 10 and will remain unused.

        - If you want to limit the plot along the horizontal axis in the
          final rendered figure, then pass the ``xmin`` and ``xmax``
          keywords to the :meth:`~sage.plot.graphics.Graphics.show` method.
          To limit the plot along the vertical axis, ``ymin`` and ``ymax``
          keywords can be provided to either this ``plot`` command or to
          the ``show`` command.

        - For the other keyword options that the ``plot`` function can
          take, refer to the method :meth:`~sage.plot.graphics.Graphics.show`.

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

    - ``linestyle`` - (default: "-") The style of the line, which is one of
       - ``"-"`` or ``"solid"``
       - ``"--"`` or ``"dashed"``
       - ``"-."`` or ``"dash dot"``
       - ``":"`` or ``"dotted"``
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

      - a dictionary ``d`` (only if a list of functions are plotted):
        The keys of the dictionary should be integers.
        The value of ``d[i]`` specifies the fill options for the i-th function
        in the list.  If ``d[i] == [j]``: Fill the area between the i-th and
        the j-th function in the list.  (But if ``d[i] == j``: Fill the area
        between the i-th function in the list and the horizontal line y = j.)

    - ``fillcolor`` - (default: 'automatic') The color of the fill.
      Either 'automatic' or a color.

    - ``fillalpha`` - (default: 0.5) How transparent the fill is.
      A number between 0 and 1.

    Note that this function does NOT simply sample equally spaced
    points between ``xmin`` and ``xmax``. Instead it computes equally spaced
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

    It is possible to use a different color for the text of each label::

        sage: p1 = plot(sin, legend_label='sin', legend_color='red')
        sage: p2 = plot(cos, legend_label='cos', legend_color='green')
        sage: p1 + p2

    Note that the independent variable may be omitted if there is no
    ambiguity::

        sage: plot(sin(1/x), (-1, 1))

    Plotting in logarithmic scale is possible for 2D plots.  There
    are two different syntaxes supported::

        sage: plot(exp, (1, 10), scale='semilogy') # log axis on vertical

    ::

        sage: plot_semilogy(exp, (1, 10)) # same thing

    ::

        sage: plot_loglog(exp, (1, 10))   # both axes are log

    ::

        sage: plot(exp, (1, 10), scale='loglog', base=2) # long time # base of log is 2

    We can also change the scale of the axes in the graphics just before
    displaying::

        sage: G = plot(exp, 1, 10) # long time
        sage: G.show(scale=('semilogy', 2)) # long time

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

    Another way to avoid getting complex numbers for negative input is to
    calculate for the positive and negate the answer::

        sage: plot(sign(x)*abs(x)^(1/3),-1,1)

    We can detect the poles of a function::

        sage: plot(gamma, (-3, 4), detect_poles = True).show(ymin = -5, ymax = 5)

    We draw the Gamma-Function with its poles highlighted::

        sage: plot(gamma, (-3, 4), detect_poles = 'show').show(ymin = -5, ymax = 5)

    The basic options for filling a plot::

        sage: p1 = plot(sin(x), -pi, pi, fill = 'axis')
        sage: p2 = plot(sin(x), -pi, pi, fill = 'min')
        sage: p3 = plot(sin(x), -pi, pi, fill = 'max')
        sage: p4 = plot(sin(x), -pi, pi, fill = 0.5)
        sage: graphics_array([[p1, p2], [p3, p4]]).show(frame=True, axes=False) # long time

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
    must use dictionary key-value pairs ``i:[j]``; key-value pairs
    like ``i:j`` will fill between the ith function and the line y=j::

        sage: def b(n): return lambda x: bessel_J(n, x) + 0.5*(n-1)
        sage: plot([b(c) for c in [1..5]], 0, 40, fill = dict([(i, [i+1]) for i in [0..3]]))
        sage: plot([b(c) for c in [1..5]], 0, 40, fill = dict([(i, i+1) for i in [0..3]])) # long time

    Extra options will get passed on to :meth:`~sage.plot.graphics.Graphics.show`,
    as long as they are valid::

        sage: plot(sin(x^2), (x, -3, 3), title='Plot of $\sin(x^2)$', axes_labels=['$x$','$y$']) # These labels will be nicely typeset
        sage: plot(sin(x^2), (x, -3, 3), title='Plot of sin(x^2)', axes_labels=['x','y']) # These will not

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

    You can even have custom tick labels along with custom positioning. ::

        sage: plot(x**2, (x,0,3), ticks=[[1,2.5],[0.5,1,2]], tick_formatter=[["$x_1$","$x_2$"],["$y_1$","$y_2$","$y_3$"]])

    You can force Type 1 fonts in your figures by providing the relevant
    option as shown below. This also requires that LaTeX, dvipng and
    Ghostscript be installed::

        sage: plot(x, typeset='type1') # optional - latex

    A example with excluded values::

        sage: plot(floor(x), (x, 1, 10), exclude = [1..10])

    We exclude all points where :class:`~sage.functions.prime_pi.PrimePi`
    makes a jump::

        sage: jumps = [n for n in [1..100] if prime_pi(n) != prime_pi(n-1)]
        sage: plot(lambda x: prime_pi(x), (x, 1, 100), exclude = jumps)

    Excluded points can also be given by an equation::

        sage: g(x) = x^2-2*x-2
        sage: plot(1/g(x), (x, -3, 4), exclude = g(x) == 0, ymin = -5, ymax = 5) # long time 

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

    We test that we can plot `f(x)=x` (see :trac:`10246`)::

        sage: f(x)=x; f
        x |--> x
        sage: plot(f,(x,-1,1))

    Check that :trac:`15030` is fixed::

        sage: plot(abs(log(x)), x)
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
    """
    Internal function which does the actual plotting.

    INPUT:

    - ``funcs`` - function or list of functions to be plotted
    - ``xrange`` - two or three tuple of [input variable], min and max
    - ``parametric`` - (default: False) a boolean for whether
      this is a parametric plot
    - ``polar`` - (default: False) a boolean for whether
      this is a polar plot
    - ``fill`` - (default: False) an input for whether
      this plot is filled
    - ``randomize`` - (default: True) a boolean for whether
      to use random plot points

    The following option is deprecated in favor of ``legend_label``:

    - ``label`` - (default: '') a string for the label

    All other usual plot options are also accepted, and a number
    are required (see the example below) which are normally passed
    through the options decorator to :func:`plot`.

    OUTPUT:

    - A ``Graphics`` object

    EXAMPLES::

    See :func:`plot` for many, many implicit examples.
    Here is an explicit one::

        sage: from sage.plot.plot import _plot
        sage: P = _plot(e^(-x^2),(-3,3),fill=True,color='red',plot_points=50,adaptive_tolerance=2,adaptive_recursion=True,exclude=None)
        sage: P.show(aspect_ratio='automatic')

    TESTS:

    Make sure that we get the right number of legend entries as the number of
    functions varies (:trac:`10514`)::

        sage: p1 = plot(1*x, legend_label='1x')
        sage: p2 = plot(2*x, legend_label='2x', color='green')
        sage: p1+p2

    ::

        sage: len(p1.matplotlib().axes[0].legend().texts)
        1
        sage: len((p1+p2).matplotlib().axes[0].legend().texts)
        2
        sage: q1 = plot([sin(x), tan(x)], legend_label='trig')
        sage: len((q1).matplotlib().axes[0].legend().texts) # used to raise AttributeError
        1
        sage: q1

    ::

    Make sure that we don't get multiple legend labels for plot segments
    (:trac:`11998`)::

        sage: p1 = plot(1/(x^2-1),(x,-2,2),legend_label="foo",detect_poles=True)
        sage: len(p1.matplotlib().axes[0].legend().texts)
        1
        sage: p1.show(ymin=-10,ymax=10) # should be one legend

    """

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
        from sage.plot.colors import rainbow
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
            if i >= 1:
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
    legend_label = options.pop('legend_label', None)
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

        G += line(data[start_index:], legend_label=legend_label, **options)
    else:
        G += line(data, legend_label=legend_label, **options)

    # Label?
    if label:
        from sage.misc.superseded import deprecation
        deprecation(4342, "Consider using legend_label instead")
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

        sage: parametric_plot( (5*cos(x), 5*sin(x), x), (x,-12, 12), plot_points=150, color="red") # long time

        sage: y=var('y')
        sage: parametric_plot( (5*cos(x), x*y, cos(x*y)), (x, -4,4), (y,-4,4)) # long time`

        sage: t=var('t')
        sage: parametric_plot( vector((sin(t), sin(2*t))), (t, 0, 2*pi), color='green') # long time
        sage: parametric_plot( vector([t, t+1, t^2]), (t, 0, 1)) # long time

    Plotting in logarithmic scale is possible with 2D plots. The keyword
    ``aspect_ratio`` will be ignored if the scale is not ``'loglog'`` or
    ``'linear'``.::

        sage: parametric_plot((x, x**2), (x, 1, 10), scale='loglog')

    We can also change the scale of the axes in the graphics just before
    displaying. In this case, the ``aspect_ratio`` must be specified as
    ``'automatic'`` if the ``scale`` is set to ``'semilogx'`` or ``'semilogy'``. For
    other values of the ``scale`` parameter, any ``aspect_ratio`` can be
    used, or the keyword need not be provided.::

        sage: p = parametric_plot((x, x**2), (x, 1, 10))
        sage: p.show(scale='semilogy', aspect_ratio='automatic')

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
        from sage.misc.superseded import deprecation
        deprecation(7008, "variable ranges to parametric_plot must be given as tuples, like (2,4) or (t,2,3)")
        args=tuple(args)
        num_ranges=1

    num_funcs = len(funcs)

    num_vars=len(sage.plot.misc.unify_arguments(funcs)[0])
    if num_vars>num_ranges:
        raise ValueError("there are more variables than variable ranges")

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = kwargs.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        kwargs['aspect_ratio'] = 'automatic'

    if num_funcs == 2 and num_ranges == 1:
        kwargs['parametric'] = True
        return plot(funcs, *args, **kwargs)
    elif (num_funcs == 3 and num_ranges <= 2):
        return sage.plot.plot3d.parametric_plot3d.parametric_plot3d(funcs, *args, **kwargs)
    else:
        raise ValueError("the number of functions and the number of variable ranges is not a supported combination for a 2d or 3d parametric plots")

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
    ``list_plot`` takes either a list of numbers, a list of tuples, a numpy
    array, or a dictionary and plots the corresponding points.

    If given a list of numbers (that is, not a list of tuples or lists),
    ``list_plot`` forms a list of tuples ``(i, x_i)`` where ``i`` goes from
    0 to ``len(data)-1`` and ``x_i`` is the ``i``-th data value, and puts
    points at those tuple values.

    ``list_plot`` will plot a list of complex numbers in the obvious
    way; any numbers for which
    :func:`CC()<sage.rings.complex_field.ComplexField>` makes sense will
    work.

    ``list_plot`` also takes a list of tuples ``(x_i, y_i)`` where ``x_i``
    and ``y_i`` are the ``i``-th values representing the ``x``- and
    ``y``-values, respectively.

    If given a dictionary, ``list_plot`` interprets the keys as
    `x`-values and the values as `y`-values.

    The ``plotjoined=True`` option tells ``list_plot`` to plot a line
    joining all the data.

    It is possible to pass empty dictionaries, lists, or tuples to
    ``list_plot``. Doing so will plot nothing (returning an empty plot).

    EXAMPLES::

        sage: list_plot([i^2 for i in range(5)]) # long time

    Here are a bunch of random red points::

        sage: r = [(random(),random()) for _ in range(20)]
        sage: list_plot(r,color='red')

    This gives all the random points joined in a purple line::

        sage: list_plot(r, plotjoined=True, color='purple')

    You can provide a numpy array.::

        sage: import numpy
        sage: list_plot(numpy.arange(10))

        sage: list_plot(numpy.array([[1,2], [2,3], [3,4]]))

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

    Plotting in logarithmic scale is possible for 2D list plots.
    There are two different syntaxes available::

        sage: yl = [2**k for k in range(20)]
        sage: list_plot(yl, scale='semilogy')  # long time  # log axis on vertical

    ::

        sage: list_plot_semilogy(yl)       # same

    .. warning::

        If ``plotjoined`` is ``False`` then the axis that is in log scale
        must have all points strictly positive. For instance, the following
        plot will show no points in the figure since the points in the
        horizontal axis starts from `(0,1)`.

        ::

            sage: list_plot(yl, scale='loglog')         # both axes are log

        Instead this will work. We drop the point `(0,1)`.::

            sage: list_plot(zip(range(1,len(yl)), yl[1:]), scale='loglog') # long time 

    We use :func:`list_plot_loglog` and plot in a different base.::

        sage: list_plot_loglog(zip(range(1,len(yl)), yl[1:]), base=2) # long time

    We can also change the scale of the axes in the graphics just before
    displaying::

        sage: G = list_plot(yl) # long time
        sage: G.show(scale=('semilogy', 2)) # long time

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
    try:
        if not data:
            return Graphics()
    except ValueError: # numpy raises ValueError if it is not empty
        pass
    if not isinstance(plotjoined, bool):
        raise TypeError("The second argument 'plotjoined' should be boolean "
                    "(True or False).  If you meant to plot two lists 'x' "
                    "and 'y' against each other, use 'list_plot(zip(x,y))'.")
    if isinstance(data, dict):
        if plotjoined:
            list_data = sorted(list(data.iteritems()))
        else:
            list_data = list(data.iteritems())
        return list_plot(list_data, plotjoined=plotjoined, **kwargs)
    try:
        from sage.rings.all import RDF
        tmp = RDF(data[0])
        data = list(enumerate(data))
    except TypeError:
        pass
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
        # if we get here, we already did "list(enumerate(data))",
        # so look at z[1] in inner list
        data = [(z.real(), z.imag()) for z in [CC(z[1]) for z in data]]
        if plotjoined:
            return line(data, **kwargs)
        else:
            return point(data, **kwargs)

#------------------------ Graphs on log scale ---------------------------#
@options(base=10)
def plot_loglog(funcs, *args, **kwds):
    """
    Plot graphics in 'loglog' scale, that is, both the horizontal and the
    vertical axes will be in logarithmic scale.

    INPUTS:

    - ``base`` -- (default: 10) the base of the logarithm. This must be
      greater than 1. The base can be also given as a list or tuple
      ``(basex, basey)``.  ``basex`` sets the base of the logarithm along the
      horizontal axis and ``basey`` sets the base along the vertical axis.

    - ``funcs`` -- any Sage object which is acceptable to the :func:`plot`.

    For all other inputs, look at the documentation of :func:`plot`.

    EXAMPLES::

        sage: plot_loglog(exp, (1,10)) # plot in loglog scale with base 10

    ::

        sage: plot_loglog(exp, (1,10), base=2.1) # long time # with base 2.1 on both axes

    ::

        sage: plot_loglog(exp, (1,10), base=(2,3))

    """
    return plot(funcs, *args, scale='loglog', **kwds)

@options(base=10)
def plot_semilogx(funcs, *args, **kwds):
    """
    Plot graphics in 'semilogx' scale, that is, the horizontal axis will be
    in logarithmic scale.

    INPUTS:

    - ``base`` -- (default: 10) the base of the logarithm. This must be
      greater than 1.

    - ``funcs`` -- any Sage object which is acceptable to the :func:`plot`.

    For all other inputs, look at the documentation of :func:`plot`.

    EXAMPLES::

        sage: plot_semilogx(exp, (1,10)) # long time # plot in semilogx scale, base 10

    ::

        sage: plot_semilogx(exp, (1,10), base=2) # with base 2

    """
    return plot(funcs, *args, scale='semilogx', **kwds)

@options(base=10)
def plot_semilogy(funcs, *args, **kwds):
    """
    Plot graphics in 'semilogy' scale, that is, the vertical axis will be
    in logarithmic scale.

    INPUTS:

    - ``base`` -- (default: 10) the base of the logarithm. This must be
      greater than 1.

    - ``funcs`` -- any Sage object which is acceptable to the :func:`plot`.

    For all other inputs, look at the documentation of :func:`plot`.

    EXAMPLES::

        sage: plot_semilogy(exp, (1,10)) # long time # plot in semilogy scale, base 10

    ::

        sage: plot_semilogy(exp, (1,10), base=2) # long time # with base 2

    """
    return plot(funcs, *args, scale='semilogy', **kwds)

@options(base=10)
def list_plot_loglog(data, plotjoined=False, **kwds):
    """
    Plot the ``data`` in 'loglog' scale, that is, both the horizontal and the
    vertical axes will be in logarithmic scale.

    INPUTS:

    - ``base`` -- (default: 10) the base of the logarithm. This must be
      greater than 1. The base can be also given as a list or tuple
      ``(basex, basey)``.  ``basex`` sets the base of the logarithm along the
      horizontal axis and ``basey`` sets the base along the vertical axis.

    For all other inputs, look at the documentation of :func:`list_plot`.


    EXAMPLES::

        sage: yl = [5**k for k in range(10)]; xl = [2**k for k in range(10)]
        sage: list_plot_loglog(zip(xl, yl)) # long time # plot in loglog scale with base 10

    ::

        sage: list_plot_loglog(zip(xl, yl), base=2.1) # long time # with base 2.1 on both axes

    ::

        sage: list_plot_loglog(zip(xl, yl), base=(2,5)) # long time

    .. warning::

        If ``plotjoined`` is ``False`` then the axis that is in log scale
        must have all points strictly positive. For instance, the following
        plot will show no points in the figure since the points in the
        horizontal axis starts from `(0,1)`.

        ::

            sage: yl = [2**k for k in range(20)]
            sage: list_plot_loglog(yl)

        Instead this will work. We drop the point `(0,1)`.::

            sage: list_plot_loglog(zip(range(1,len(yl)), yl[1:]))

    """
    return list_plot(data, plotjoined=plotjoined, scale='loglog', **kwds)

@options(base=10)
def list_plot_semilogx(data, plotjoined=False, **kwds):
    """
    Plot ``data`` in 'semilogx' scale, that is, the horizontal axis will be
    in logarithmic scale.

    INPUTS:

    - ``base`` -- (default: 10) the base of the logarithm. This must be
      greater than 1.

    For all other inputs, look at the documentation of :func:`list_plot`.

    EXAMPLES::

        sage: yl = [2**k for k in range(12)]
        sage: list_plot_semilogx(zip(yl,yl))

    .. warning::

        If ``plotjoined`` is ``False`` then the horizontal axis must have all
        points strictly positive. Otherwise the plot will come up empty.
        For instance the following plot contains a point at `(0,1)`.

        ::

            sage: yl = [2**k for k in range(12)]
            sage: list_plot_semilogx(yl) # plot is empty because of `(0,1)`

        We remove `(0,1)` to fix this.::

            sage: list_plot_semilogx(zip(range(1, len(yl)), yl[1:]))

    ::

        sage: list_plot_semilogx([(1,2),(3,4),(3,-1),(25,3)], base=2) # with base 2

    """
    return list_plot(data, plotjoined=plotjoined, scale='semilogx', **kwds)

@options(base=10)
def list_plot_semilogy(data, plotjoined=False, **kwds):
    """
    Plot ``data`` in 'semilogy' scale, that is, the vertical axis will be
    in logarithmic scale.

    INPUTS:

    - ``base`` -- (default: 10) the base of the logarithm. This must be
      greater than 1.

    For all other inputs, look at the documentation of :func:`list_plot`.

    EXAMPLES::

        sage: yl = [2**k for k in range(12)]
        sage: list_plot_semilogy(yl) # plot in semilogy scale, base 10

    .. warning::

        If ``plotjoined`` is ``False`` then the vertical axis must have all
        points strictly positive. Otherwise the plot will come up empty.
        For instance the following plot contains a point at `(1,0)`.

        ::

            sage: xl = [2**k for k in range(12)]; yl = range(len(xl))
            sage: list_plot_semilogy(zip(xl,yl)) # plot empty due to (1,0)

        We remove `(1,0)` to fix this.::

            sage: list_plot_semilogy(zip(xl[1:],yl[1:]))


    ::

        sage: list_plot_semilogy([2, 4, 6, 8, 16, 31], base=2) # with base 2

    """
    return list_plot(data, plotjoined=plotjoined, scale='semilogy', **kwds)

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


def reshape(v, n, m):
    """
    Helper function for creating graphics arrays.

    The input array is flattened and turned into an `n\times m`
    array, with blank graphics object padded at the end, if
    necessary.

    INPUT:

    -  ``v`` - a list of lists or tuples

    -  ``n, m`` - integers

    OUTPUT:

    A list of lists of graphics objects

    EXAMPLES::

        sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in range(10)]
        sage: graphics_array(L,3,4)  # long time (up to 4s on sage.math, 2012)

    ::

        sage: M = [[plot(sin(k*x),(x,-pi,pi)) for k in range(3)],[plot(cos(j*x),(x,-pi,pi)) for j in [3..5]]]
        sage: graphics_array(M,6,1)  # long time (up to 4s on sage.math, 2012)

    TESTS::

        sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in [1..3]]
        sage: graphics_array(L,0,-1) # indirect doctest
        Traceback (most recent call last):
        ...
        AssertionError: array sizes must be positive
    """
    assert n>0 and m>0, 'array sizes must be positive'
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
        sage: p1 = plot(f,(-2*pi,2*pi),color=hue(0.5)) # long time
        sage: p2 = plot(g,(-2*pi,2*pi),color=hue(0.9)) # long time
        sage: p3 = parametric_plot((f,g),(0,2*pi),color=hue(0.6)) # long time
        sage: p4 = parametric_plot((f,h),(0,2*pi),color=hue(1.0)) # long time

    Now make a graphics array out of the plots::

        sage: graphics_array(((p1,p2),(p3,p4))) # long time

    One can also name the array, and then use :meth:`~sage.plot.graphics.GraphicsArray.show`
    or :meth:`~sage.plot.graphics.GraphicsArray.save`::

        sage: ga = graphics_array(((p1,p2),(p3,p4))) # long time
        sage: ga.show() # long time

    Here we give only one row::

        sage: p1 = plot(sin,(-4,4))
        sage: p2 = plot(cos,(-4,4))
        sage: g = graphics_array([p1, p2]); print g
        Graphics Array of size 1 x 2
        sage: g.show()

    It is possible to use ``figsize`` to change the size of the plot
    as a whole::

        sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in [1..3]]
        sage: G = graphics_array(L)
        sage: G.show(figsize=[5,3])  # smallish and compact

    ::

        sage: G.show(figsize=[10,20])  # bigger and tall and thin; long time (2s on sage.math, 2012)

    ::

        sage: G.show(figsize=8)  # figure as a whole is a square
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
        See http://trac.sagemath.org/7008 for details.
        (theta, [2.0, 2.75, 3.5, 4.25, 5.0])
        sage: var_and_list_of_values((2, 5),  5)
        (None, [2.0, 2.75, 3.5, 4.25, 5.0])
        sage: var_and_list_of_values((var('theta'), 2, 5),  2)
        (theta, [2.0, 5.0])
        sage: var_and_list_of_values((2, 5),  2)
        (None, [2.0, 5.0])
    """
    from sage.misc.superseded import deprecation
    deprecation(7008, "var_and_list_of_values is deprecated.  Please use sage.plot.misc.setup_for_eval_on_grid; note that that function has slightly different calling and return conventions which make it more generally applicable")
    plot_points = int(plot_points)
    if plot_points < 2:
        raise ValueError("plot_points must be greater than 1")
    if not isinstance(v, (tuple, list)):
        raise TypeError("v must be a tuple or list")
    if len(v) == 3:
        var = v[0]
        a, b = v[1], v[2]
    elif len(v) == 2:
        var = None
        a, b = v
    else:
        raise ValueError("parametric value range must be a list or tuple of length 2 or 3.")

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
    ``sage.plot.misc.setup_for_eval_on_grid`` instead.  Please note that
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
        See http://trac.sagemath.org/7008 for details.
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
    from sage.misc.superseded import deprecation
    deprecation(7008, "sage.plot.plot.setup_for_eval_on_grid is deprecated.  Please use sage.plot.misc.setup_for_eval_on_grid; note that that function has slightly different calling and return conventions which make it more generally applicable")

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
    The adaptive refinement algorithm for plotting a function ``f``. See
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
       for more information.  See the documentation for :func:`plot` for more
       information on how the adaptive refinement algorithm works.

    OUTPUT:


    -  ``list`` - a list of points to insert between ``p1`` and
       ``p2`` to get a better linear approximation between them


    TESTS::

        sage: from sage.plot.plot import adaptive_refinement
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01, adaptive_recursion=0)
        []
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01)
        [(0.125*pi, 0.3826834323650898), (0.1875*pi, 0.5555702330196022), (0.25*pi, 0.7071067811865475), (0.3125*pi, 0.8314696123025452), (0.375*pi, 0.9238795325112867), (0.4375*pi, 0.9807852804032304), (0.5*pi, 1.0), (0.5625*pi, 0.9807852804032304), (0.625*pi, 0.9238795325112867), (0.6875*pi, 0.8314696123025455), (0.75*pi, 0.7071067811865476), (0.8125*pi, 0.5555702330196022), (0.875*pi, 0.3826834323650899)]

    This shows that lowering ``adaptive_tolerance`` and raising
    ``adaptive_recursion`` both increase the number of subdivision
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

    except (ZeroDivisionError, TypeError, ValueError, OverflowError) as msg:
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

        except (ArithmeticError, TypeError, ValueError) as msg:
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
                    except (ArithmeticError, TypeError, ValueError) as msg:
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
                    except (ArithmeticError, TypeError, ValueError) as msg:
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

