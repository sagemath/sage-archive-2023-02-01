r"""
2D Plotting

\sage provides extensive 2D plotting functionality.  The underlying
rendering is done using the matplotlib Python library.

The following graphics primitives are supported:
\begin{itemize}
    \item arrow  -- an arrow from a min point to a max point.
    \item circle -- a circle with given radius
    \item disk   -- a filled disk
    \item line   -- a line determined by a sequence of points (this need not be straight!)
    \item point  -- a point
    \item text   -- some text
    \item polygon -- a filled polygon
\end{itemize}

The following plotting functions are supported:
\begin{itemize}
    \item plot   -- plot of a function or other \sage object (e.g., elliptic curve).
    \item parametric_plot
    \item polar_plot
    \item list_plot
    \item bar_chart
    \item contour_plot
    \item plot_vector_field
    \item matrix_plot
    \item graphics_array
\end{itemize}

The following miscellaneous Graphics functions are included:
\begin{itemize}
    \item Graphics
    \item is_Graphics
    \item rgbcolor
    \item hue
\end{itemize}

Type \kbd{?} after each primitive in \sage for help and examples.

EXAMPLES:
We construct a plot involving several graphics objects:

    sage: G = plot(cos, -5, 5, thickness=5, rgbcolor=(0.5,1,0.5))
    sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))
    sage: P   # show it

We draw a circle and a curve:
    sage: circle((1,1), 1) + plot(x^2, (0,5))

Notice that the above circle is not round, because the aspect ratio of the
coordinate system is not 1:1.    The \code{aspect_ratio} option to show
allows us to fix this:
    sage: show(circle((1,1), 1) + plot(x^2, (0,5)), aspect_ratio=1)

With an aspect ratio of 2 the circle is squashed half way down (it looks twice
as wide as it does tall):
    sage: show(circle((1,1), 1) + plot(x^2, (0,5)), aspect_ratio=2)

Use figsize to set the actual aspect ratio of the rendered image
(i.e., of the frame).  For example, this image is twice as many pixels
wide as it is tall:
    sage: show(circle((1,1), 1) + plot(x^2, (0,5)), figsize=[8,4])

Next we construct the reflection of the above polygon about the
$y$-axis by iterating over the list of first-coordinates of the first
graphic element of $P$ (which is the actual Polygon; note that $P$ is
a Graphics object, which consists of a single polygon):

    sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
    sage: Q   # show it

We combine together different graphics objects using ``+'':

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

We can put text in a graph:

    sage: L = [[cos(pi*i/100)^3,sin(pi*i/100)] for i in range(200)]
    sage: p = line(L, rgbcolor=(1/4,1/8,3/4))
    sage: t = text('A Bulb', (1.5, 0.25))
    sage: x = text('x axis', (1.5,-0.2))
    sage: y = text('y axis', (0.4,0.9))
    sage: g = p+t+x+y
    sage: g.show(xmin=-1.5, xmax=2, ymin=-1, ymax=1)

We plot the Riemann zeta function along the critical line and
see the first few zeros:

    sage: i = CDF.0      # define i this way for maximum speed.
    sage: p1 = plot(lambda t: arg(zeta(0.5+t*i)), 1,27,rgbcolor=(0.8,0,0))
    sage: p2 = plot(lambda t: abs(zeta(0.5+t*i)), 1,27,rgbcolor=hue(0.7))
    sage: print p1 + p2
    Graphics object consisting of 2 graphics primitives
    sage: p1 + p2    # display it

Many concentric circles shrinking toward the origin:
    sage: show(sum(circle((i,0), i, hue=sin(i/10)) for i in [10,9.9,..,0]), aspect_ratio=1)

Here is a pretty graph:
    sage: g = Graphics()
    sage: for i in range(60):
    ...    p = polygon([(i*cos(i),i*sin(i)), (0,i), (i,0)],\
    ...                rgbcolor=hue(i/40+0.4), alpha=0.2)
    ...    g = g + p
    ...
    sage: g.show(dpi=200, axes=False)

Another graph:
    sage: x = var('x')
    sage: P = plot(sin(x)/x, -4,4, rgbcolor=(0,0,1)) + \
    ...       plot(x*cos(x), -4,4, rgbcolor=(1,0,0)) + \
    ...       plot(tan(x),-4,4,rgbcolor=(0,1,0))
    ...
    sage: P.show(ymin=-pi,ymax=pi)

PYX EXAMPLES:
These are some examples of plots similar to some of the plots in the
PyX (\url{http://pyx.sourceforge.net}) documentation:

Symbolline:
    sage: y(x) = x*sin(x^2)
    sage: v = [(x, y(x)) for x in [-3,-2.95,..,3]]
    sage: show(points(v, rgbcolor=(0.2,0.6, 0.1), pointsize=30) + plot(spline(v), -3.1, 3))

Cycliclink:
    sage: x = var('x')
    sage: g1 = plot(cos(20*x)*exp(-2*x), 0, 1)
    sage: g2 = plot(2*exp(-30*x) - exp(-3*x), 0, 1)
    sage: show(graphics_array([g1, g2], 2, 1), xmin=0)

Pi Axis:
In the PyX manual, the point of this example is to show labeling the
X-axis using rational multiples of Pi.  Sage currently has no support
for controlling how the ticks on the x and y axes are labeled, so
this is really a bad example:

    sage: g1 = plot(sin(x), 0, 2*pi)
    sage: g2 = plot(cos(x), 0, 2*pi, linestyle = "--")
    sage: g1 + g2    # show their sum

An illustration of integration:
    sage: f = (x-3)*(x-5)*(x-7)+40
    sage: P = line([(2,0),(2,f(2))], rgbcolor=(0,0,0))
    sage: P += line([(8,0),(8,f(8))], rgbcolor=(0,0,0))
    sage: P += polygon([(2,0),(2,f(2))] + [(x, f(x)) for x in [2,2.1,..,8]] + [(8,0),(2,0)],  rgbcolor=(0.8,0.8,0.8))
    sage: P += text("$\\int_{a}^b f(x) dx$", (5, 20), fontsize=16, rgbcolor=(0,0,0))
    sage: P += plot(f, 1, 8.5, thickness=3)
    sage: P    # show the result

NUMERICAL PLOTTING:

\sage also provides 2D plotting with an interface that is a likely very
familiar to people doing numerical computation.  For example,

    sage: from pylab import *
    sage: t = arange(0.0, 2.0, 0.01)
    sage: s = sin(2*pi*t)
    sage: P = plot(t, s, linewidth=1.0)
    sage: xl = xlabel('time (s)')
    sage: yl = ylabel('voltage (mV)')
    sage: t = title('About as simple as it gets, folks')
    sage: grid(True)
    sage: savefig('sage.png')

Since the above overwrites many Sage plotting functions, we
reset the state of Sage, so that the examples below work!
    sage: reset()

See \url{http://matplotlib.sourceforge.net} for complete documentation
about how to use Matplotlib.

TESTS:
We test dumping and loading a plot.
    sage: p = plot(sin(x), (x, 0,2*pi))
    sage: Q = loads(dumps(p))



AUTHORS:
    -- Alex Clemesha and William Stein (2006-04-10): initial version
    -- David Joyner: examples
    -- Alex Clemesha (2006-05-04) major update
    -- Willaim Stein (2006-05-29): fine tuning, bug fixes, better server integration
    -- William Stein (2006-07-01): misc polish
    -- Alex Clemesha (2006-09-29): added contour_plot, frame axes, misc polishing
    -- Robert Miller (2006-10-30): tuning, NetworkX primitive
    -- Alex Clemesha (2006-11-25): added plot_vector_field, matrix_plot, arrow,
                                   bar_chart, Axes class usage (see axes.py)
    -- Bobby Moretti and William Stein (2008-01): Change plot to specify ranges
                                   using the (varname, min, max) notation.
    -- William Stein (2008-01-19): raised the documentation coverage
                                   from a miserable 12 percent to a
                                   'wopping' 35 percent, and fixed and
                                   clarified numerous small issues.
"""

############################################################################
#  Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com> and William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

import types

from sage.structure.sage_object import SageObject

## IMPORTANT: Do *not* import matplotlib at module scope.  It takes a
## surprisingliy long time to initialize itself.  It's better if it is
## imported in functions, so it only gets started if it is actually
## going to be used.

DEFAULT_FIGSIZE=(6, 3.70820393249937)
DEFAULT_DPI = 100
EMBEDDED_MODE = False
DOCTEST_MODE = False
SHOW_DEFAULT = True

def show_default(default=None):
    r"""
    Set the default for showing plots using any plot commands.
    If called with no arguments, returns the current default.

    If this is \code{True} (the default) then any plot object when displayed
    will be displayed as an actual plot instead of text, i.e., the
    show command is not needed.

    EXAMPLES:
    The default starts out as \code{True}:
        sage: show_default()
        True

    We set it to \code{False}.
        sage: show_default(False)

    We see that it is \code{False}.
        sage: show_default()
        False

    Now plot commands will not display their plots by default.

    Turn back on default display.
        sage: show_default(True)
    """
    global SHOW_DEFAULT
    if default is None:
        return SHOW_DEFAULT
    SHOW_DEFAULT = bool(default)

do_verify = True

from sage.misc.randstate import current_randstate #for plot adaptive refinement
import os #for viewing and writing images
from colorsys import hsv_to_rgb #for the hue function
from math import sin, cos, modf, pi #for hue and polar_plot
from sage.structure.sage_object import SageObject

from sage.ext.fast_eval import fast_float, fast_float_constant, is_fast_float

import sage.misc.misc

from misc import rgbcolor, Color, options, rename_keyword, to_mpl_color

import operator

############### WARNING ###
# Try not to import any matplotlib stuff here -- matplotlib is
# slow to import.  (I did benchmarking and found that by not
# importing here, and instead importing when needed below, that
# Sage startup times are much improved.)  - William
###############

#Sage 2D Graphics Axes class:
from axes import Axes
from axes import GridLines

def is_Graphics(x):
    """
    Return True if $x$ is a Graphics object.

    EXAMPLES:
        sage: from sage.plot.plot import is_Graphics
        sage: is_Graphics(1)
        False
        sage: is_Graphics(disk((0.0, 0.0), 1, (0, pi/2)))
        True
    """
    return isinstance(x, Graphics)

class Graphics(SageObject):
    """
    The Graphics object is an empty list of graphics objects
    It is useful to use this object when intializing a
    for loop where different graphics object will be added
    to the empty object.

    EXAMPLES:
        sage: G = Graphics(); print G
        Graphics object consisting of 0 graphics primitives
        sage: c = circle((1,1), 1)
        sage: G+=c; print G
        Graphics object consisting of 1 graphics primitive

    Here we make a graphic of embedded isosceles triangles,
    coloring each one with a different color as we go:

        sage: h=10; c=0.4; p=0.1;
        sage: G = Graphics()
        sage: for x in srange(1,h+1):
        ...        l = [[0,x*sqrt(3)],[-x/2,-x*sqrt(3)/2],[x/2,-x*sqrt(3)/2],[0,x*sqrt(3)]]
        ...        G+=line(l,rgbcolor=hue(c + p*(x/h)))
        sage: G.show(figsize=[5,5])

    """

    def __init__(self):
        """
        Create a new empty Graphics objects with all the defaults.

        EXAMPLES:
            sage: G = Graphics()
        """
        self.__aspect_ratio = None
        self.__fontsize = 10
        self.__show_axes = True
        self.__axes_color = (0, 0, 0)
        self.__axes_label_color = (0, 0, 0)
        self.__tick_label_color = (0, 0, 0)
        self.__axes_width = 0.8
        self.__objects = []
        self.__axes_range = {}

    def set_aspect_ratio(self, ratio):
        """
        Set the aspect ratio.

        INPUT:
            ratio  -- a positive real number

        EXAMPLES:
        We create a plot of a circle, and it doesn't look quite round:
            sage: P = circle((1,1), 1); P

        So we set the aspect ratio and now it is round:
            sage: P.set_aspect_ratio(1)
            sage: P

        Note that the aspect ratio is inherited upon addition (which takes the
        max of aspect ratios of objects whose aspect ratio has been set):
            sage: P + circle((0,0), 0.5)           # still square

        In the following example, both plots produce a circle that looks twice
        as wide as tall:
            sage: Q = circle((0,0), 0.5); Q.set_aspect_ratio(2)
            sage: P + Q
            sage: Q + P
        """
        ratio = float(ratio)
        if ratio <= 0:
            raise ValueError, "the aspect ratio must be positive"
        self.__aspect_ratio = ratio

    def aspect_ratio(self):
        """
        Get the current aspect ratio.

        OUTPUT:
            either None if the aspect ratio hasn't been set or a positive float

        EXAMPLES:
            sage: P = circle((1,1), 1)
            sage: P.aspect_ratio() is None
            True
            sage: P.set_aspect_ratio(2)
            sage: P.aspect_ratio()
            2.0
        """
        return self.__aspect_ratio

    def axes_range(self, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Set the ranges of the $x$ and $y$ axes.

        INPUT:
            xmin, xmax, ymin, ymax -- floats

        EXAMPLES:
            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.axes_range(-1, 20, 0, 2)
            sage: d = L.get_axes_range()
            sage: d['xmin'], d['xmax'], d['ymin'], d['ymax']
            (-1.0, 20.0, 0.0, 2.0)
        """
        l = locals()
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            if l[name] is not None:
                self.__axes_range[name] = float(l[name])

    def get_axes_range(self):
        """
        EXAMPLES:
            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.get_axes_range()
            {}
            sage: L.axes_range(xmin=-1)
            sage: L.get_axes_range()
            {'xmin': -1.0}
        """
        return self.__axes_range

    def fontsize(self, s=None):
        """
        Set the font size of axes labels and tick marks.

        INPUT:
            s -- integer, a font size in points.

        If called with no input, return the current fontsize.

        EXAMPLES:
            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.fontsize()
            10
            sage: L.fontsize(20)
            sage: L.fontsize()
            20

        All the numbers on the axes will be very large in this plot:
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
        Set whether or not the $x$ and $y$ axes are shown by default.

        INPUT:
            show -- bool

        If called with no input, return the current axes setting.

        EXAMPLES:
            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])

        By default the axes are displayed.
            sage: L.axes()
            True

        But we turn them off, and verify that they are off
            sage: L.axes(False)
            sage: L.axes()
            False

        Displaying L now shows a triangle but no axes.
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
            c -- an rgb color 3-tuple, where each tuple entry is a
                 float between 0 and 1

        EXAMPLES:
        We create a line, which has like everything a default axes color of black.
            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.axes_color()
            (0, 0, 0)

        We change the axes color to red and verify the change.
            sage: L.axes_color((1,0,0))
            sage: L.axes_color()
            (1.0, 0.0, 0.0)

        When we display the plot, we'll see a blue triangle and bright red axes.
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
            l -- (default: None) a list of two strings or None

        OUTPUT:
            a 2-tuple of strings

        If l is None, returns the current \code{axes_labels}, which is
        itself by default None.  The default
        labels are both empty.

        EXAMPLES:
        We create a plot and put x and y axes labels on it.
            sage: p = plot(sin(x), (x, 0, 10))
            sage: p.axes_labels(['x','y'])
            sage: p.axes_labels()
            ('x', 'y')

        Now when you plot p, you see x and y axes labels:
            sage: p
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

        The axes labels are placed at the edge of the x and y axes,
        and are not on by default (use the \code{axes_labels} command
        to set them; see the example below).  This function just changes
        their color.

        INPUT:
            c -- an rgb 3-tuple of numbers between 0 and 1

        If called with no input, return the current axes_label_color setting.

        EXAMPLES:
        We create a plot, which by default has axes label color black.
            sage: p = plot(sin, (-1,1))
            sage: p.axes_label_color()
            (0, 0, 0)

        We change the labels to be red, and confirm this:
            sage: p.axes_label_color((1,0,0))
            sage: p.axes_label_color()
            (1.0, 0.0, 0.0)

        We set labels, since otherwise we won't see anything.
            sage: p.axes_labels(['$x$ axis', '$y$ axis'])

        In the plot below, notice that the labels are red:
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
        Set the axes width.  Use this to draw a plot with really fat
        or really thin axes.

        INPUT:
            w -- a float

        If called with no input, return the current \code{axes_width} setting.

        EXAMPLE:
        We create a plot, see the default axes width (with funny Python float rounding),
        then reset the width to 10 (very fat).
            sage: p = plot(cos, (-3,3))
            sage: p.axes_width()
            0.80000000000000004
            sage: p.axes_width(10)
            sage: p.axes_width()
            10.0

        Finally we plot the result, which is a graph with very fat axes.
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
            c -- an rgb 3-tuple of numbers between 0 and 1

        If called with no input, return the current tick_label_color setting.

        EXAMPLES:
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

        If the \code{show_default} function has been called with True
        (the default), then you'll see this graphics object displayed.
        Otherwise you'll see a text representation of it.

        EXAMPLES:
        We create a plot and call \code{_repr_} on it, which causes it
        to be displayed as a plot:
            sage: P = plot(cos, (-1,1))
            sage: P._repr_()
            ''

        Just doing this also displays the plot:
            sage: P

        Note that printing P with the \code{print} statement does not display the plot:
            sage: print P
            Graphics object consisting of 1 graphics primitive

        Now we turn off showing plots by default:
            sage: show_default(False)

        Now we just get a string.  To show P you would have to do \code{show(P)}.
            sage: P._repr_()
            'Graphics object consisting of 1 graphics primitive'
            sage: P
            Graphics object consisting of 1 graphics primitive

        Finally, we turn \code{show_default} back on:
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

        EXAMPLES:
            sage: S = circle((0,0), 2); S.__str__()
            'Graphics object consisting of 1 graphics primitive'
            sage: print S
            Graphics object consisting of 1 graphics primitive

        WARNING: \code{__str__} is not called when printing lists of graphics
        objects, which can be confusing, since they will all pop up.  One
        workaround is to call \code{show_default}:

        For example, below when we do \code{print v} two plots are displayed:
            sage: v = [circle((0,0), 2), circle((2,3), 1)]
            sage: print v
            [, ]

        However, if we call \code{show_default} then we see the text representations
        of the graphics:
            sage: show_default(False)
            sage: print v
            [Graphics object consisting of 1 graphics primitive, Graphics object consisting of 1 graphics primitive]
            sage: v
            [Graphics object consisting of 1 graphics primitive,
             Graphics object consisting of 1 graphics primitive]

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

        EXAMPLE:
            sage: G = circle((1,1),2) + circle((2,2),5); print G
            Graphics object consisting of 2 graphics primitives
            sage: G[1]
            Circle defined by (2.0,2.0) with r=5.0
        """
        return self.__objects[i]

    def __len__(self):
        """
        If G is of type Graphics, then len(G)
        gives the number of distinct graphics
        primitives making up that object.

        EXAMPLES:
            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
        """
        return len(self.__objects)

    def __delitem__(self, i):
        """
        If G is of type Graphics, then del(G[i])
        removes the ith distinct graphic
        primitive making up that object.

        EXAMPLES:
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
        You can replace a GraphicPrimitive (point, line, circle, etc...)
        in a Graphics object G with any other GraphicPrimitive

        EXAMPLES:
            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print G
            Graphics object consisting of 3 graphics primitives

            sage: p = polygon([[1,3],[2,-2],[1,1],[1,3]]); print p
            Graphics object consisting of 1 graphics primitive

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

        This only works when other is a Python int equal to 0. In all
        other cases a TypeError is raised.  The main reason for this
        function is to make suming a list of graphics objects easier.

        EXAMPLES:
            sage: S = circle((0,0), 2)
            sage: print int(0) + S
            Graphics object consisting of 1 graphics primitive
            sage: print S + int(0)
            Graphics object consisting of 1 graphics primitive

        The following would fail were it not for this function:
            sage: v = [circle((0,0), 2), circle((2,3), 1)]
            sage: print sum(v)
            Graphics object consisting of 2 graphics primitives
        """
        if isinstance(other, (int, long)) and other == 0:
            return self
        raise TypeError

    def __add__(self, other):
        """
        If you have any Graphics object G1, you can always add any
        other amount of Graphics objects G2,G3,...  to form a new
        Graphics object: G4 = G1 + G2 + G3.

        The xmin, xmax, ymin, and ymax properties of the graphics objects
        are expanded to include all objects in both scenes.  If the aspect
        ratio property of either or both objects are set, then the larger
        aspect ratio is chosen.

        EXAMPLES:
            sage: g1 = plot(abs(sqrt(x^3-1)), (x,1,5))
            sage: g2 = plot(-abs(sqrt(x^3-1)), (x,1,5), rgbcolor=(1,0,0))
            sage: g1 + g2  # displays the plot
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
        g.__aspect_ratio = max(self.__aspect_ratio, other.__aspect_ratio)
        return g

    def add_primitive(self, primitive):
        """
        Adds a primitive to this graphics object.
        """
        self.__objects.append(primitive)

    def plot(self, *args, **kwds):
        """
        Draw a 2d plot this graphics object, which just returns this
        object since this is already a 2d graphics object.

        EXAMPLES:
            sage: S = circle((0,0), 2)
            sage: S.plot() is S
            True
        """
        return self

    def plot3d(self, z=0, **kwds):
        """
        Returns an embedding of this 2D plot into the xy-plane of 3D space, as
        a 3D plot object. An optional parameter z can be given to specify the
        z-coordinate.

        EXAMPLES:
            sage: sum([plot(z*sin(x), 0, 10).plot3d(z) for z in range(6)]) #long
        """
        from sage.plot.plot3d.base import Graphics3dGroup
        g = Graphics3dGroup([g.plot3d(**kwds) for g in self.__objects])
        if z:
            g = g.translate(0,0,z)
        return g

    def show(self, xmin=None, xmax=None, ymin=None, ymax=None,
             figsize=DEFAULT_FIGSIZE, filename=None,
             dpi=DEFAULT_DPI, axes=None, axes_labels=None,frame=False,
             fontsize=None, aspect_ratio=None,
             gridlines=None, gridlinesstyle=None,
             vgridlinesstyle=None, hgridlinesstyle=None):
        """
        Show this graphics image with the default image viewer.

        OPTIONAL INPUT:
            filename     -- (default: None) string
            dpi          -- dots per inch
            figsize      -- [width, height]
            aspect_ratio -- the perceived width divided by the
                            perceived height.  If the aspect ratio is
                            set to 1, circles will look round.  If it
                            is set to 2 they will look twice as wide
                            as they are tall.  This is the
                            aspect_ratio of the image, not of the
                            frame that contains it.  If you want to
                            set the aspect ratio of the frame, use
                            figsize.
            axes         -- (default: True)
            axes_labels  -- (default: None) list (or tuple) of two strings;
                            the first is used as the label for the horizontal
                            axis, and the second for the vertical axis.
            fontsize     -- (default: current setting -- 10) positive
                            integer; used for axes labels; if you make
                            this very large, you may have to increase
                            figsize to see all labels.
            frame        -- (default: False) draw a frame around the image
            gridlines    -- (default: None) can be any of the following:
                            1. None, False: do not add grid lines.
                            2. True, "automatic", "major": add grid lines
                               at major ticks of the axes.
                            3. "minor": add grid at major and minor ticks.
                            4. [xlist,ylist]: a tuple or list containing
                               two elements, where xlist (or ylist) can be
                               any of the following.
                               4a. None, False: don't add horizontal (or
                                   vertical) lines.
                               4b. True, "automatic", "major": add
                                   horizontal (or vertical) grid lines at
                                   the major ticks of the axes.
                               4c. "minor": add horizontal (or vertical)
                                   grid lines at major and minor ticks of
                                   axes.
                               4d. an iterable yielding numbers n or pairs
                                   (n,opts), where n is the coordinate of
                                   the line and opt is a dictionary of
                                   MATPLOTLIB options for rendering the
                                   line.
            gridlinesstyle,
            hgridlinesstyle,
            vgridlinesstyle
                         -- (default: None) a dictionary of MATPLOTLIB
                            options for the rendering of the grid lines,
                            the horizontal grid lines or the vertical grid
                            lines, respectively.

        EXAMPLES:
            sage: c = circle((1,1), 1, rgbcolor=(1,0,0))
            sage: c.show(xmin=-1, xmax=3, ymin=-1, ymax=3)

        To correct the apect ratio of certain graphics, it is necessary
        to show with a `\code{figsize}' of square dimensions.

            sage: c.show(figsize=[5,5], xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can turn off the drawing of the axes:

            sage: show(plot(sin,-4,4), axes=False)

        You can also label the axes:

            sage: show(plot(sin,-4,4), axes_labels=('x','y'))

        You can turn on the drawing of a frame around the plots:

            sage: show(plot(sin,-4,4), frame=True)

        Add grid lines at the major ticks of the axes.
            sage: c = circle((0,0), 1)
            sage: c.show(gridlines=True)
            sage: c.show(gridlines="automatic")
            sage: c.show(gridlines="major")

        Add grid lines at the major and minor ticks of the axes.
            sage: u,v = var('u v')
            sage: f = exp(-(u^2+v^2))
            sage: p = plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2))
            sage: p.show(gridlines="minor")

        Add only horizontal or vertical grid lines.
            sage: p = plot(sin,-10,20)
            sage: p.show(gridlines=[None, "automatic"])
            sage: p.show(gridlines=["minor", False])

        Add grid lines at specific positions (using lists/tuples).
            sage: x, y = var('x, y')
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3)-4*(x^2+y^2-2*x)^2, \
            ...             (-2,2), (-2,2), plot_points=1000)
            sage: p.show(gridlines=[[1,0],[-1,0,1]])

        Add grid lines at specific positions (using iterators).
            sage: def maple_leaf(t):
            ...     return (100/(100+(t-pi/2)^8))*(2-sin(7*t)-cos(30*t)/2)
            sage: p = polar_plot(maple_leaf, -pi/4, 3*pi/2, rgbcolor="red",plot_points=1000) #long
            sage: p.show(gridlines=( [-3,-2.75,..,3], xrange(-1,5,2) )) #long

        Add grid lines at specific positions (using functions).
            sage: y = x^5 + 4*x^4 - 10*x^3 - 40*x^2 + 9*x + 36
            sage: p = plot(y, -4.1, 1.1)
            sage: xlines = lambda a,b: [z for z,m in y.roots()]
            sage: p.show(gridlines=[xlines, [0]], frame=True, axes=False)

        Change the style of all the grid lines.
            sage: b = bar_chart([-3,5,-6,11], rgbcolor=(1,0,0))
            sage: b.show(gridlines=([-1,-0.5,..,4],True), \
            ...     gridlinesstyle=dict(color="blue", linestyle=":"))

        Change the style of the horizontal or vertical grid lines separately.
            sage: p = polar_plot(2 + 2*cos(x), 0, 2*pi, rgbcolor=hue(0.3))
            sage: p.show(gridlines=True, \
            ...     hgridlinesstyle=dict(color="orange", linewidth=1.0), \
            ...     vgridlinesstyle=dict(color="blue", linestyle=":"))

        Change the style of each grid line individually.
            sage: x, y = var('x, y')
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3)-4*(x^2+y^2-2*x)^2, \
            ...             (-2,2), (-2,2), plot_points=1000)
            sage: p.show(gridlines=(
            ...    [
            ...     (1,{"color":"red","linestyle":":"}),
            ...     (0,{"color":"blue","linestyle":"--"})
            ...    ],
            ...    [
            ...     (-1,{"rgbcolor":"red","linestyle":":"}),
            ...     (0,{"color":"blue","linestyle":"--"}),
            ...     (1,{"rgbcolor":"red","linestyle":":"}),
            ...    ]
            ...    ),
            ...    gridlinesstyle=dict(marker='x',rgbcolor="black"))

        Grid lines can be added to contour plots.
            sage: f = sin(x^2 + y^2)*cos(x)*sin(y)
            sage: c = contour_plot(f, (-4, 4), (-4, 4), plot_points=100)
            sage: c.show(gridlines=True, gridlinesstyle={'linestyle':':','linewidth':1, 'rgbcolor':'red'})

        Grid lines can be added to matrix plots.
            sage: M = MatrixSpace(QQ,10).random_element()
            sage: matrix_plot(M).show(gridlines=True)
        """
        if DOCTEST_MODE:
            self.save(sage.misc.misc.SAGE_TMP + '/test.png',
                      xmin, xmax, ymin, ymax, figsize,
                      dpi=dpi, axes=axes, axes_labels=axes_labels,frame=frame,
                      aspect_ratio=aspect_ratio, gridlines=gridlines,
                      gridlinesstyle=gridlinesstyle,
                      vgridlinesstyle=vgridlinesstyle,
                      hgridlinesstyle=hgridlinesstyle)
            return
        if EMBEDDED_MODE:
            self.save(filename, xmin, xmax, ymin, ymax, figsize,
                      dpi=dpi, axes=axes, axes_labels=axes_labels,frame=frame,
                      aspect_ratio=aspect_ratio, gridlines=gridlines,
                      gridlinesstyle=gridlinesstyle,
                      vgridlinesstyle=vgridlinesstyle,
                      hgridlinesstyle=hgridlinesstyle)
            return
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self.save(filename, xmin, xmax, ymin, ymax, figsize, dpi=dpi, axes=axes,
                  axes_labels=axes_labels,
                  frame=frame, fontsize=fontsize,
                  aspect_ratio=aspect_ratio,
                  gridlines=gridlines,
                  gridlinesstyle=gridlinesstyle,
                  vgridlinesstyle=vgridlinesstyle,
                  hgridlinesstyle=hgridlinesstyle)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(sage.misc.viewer.browser(), filename))


    def get_minmax_data(self):
        objects = self.__objects
        if objects:
            minmax_data = [o.get_minmax_data() for o in objects]
            xmin = min(d['xmin'] for d in minmax_data)
            xmax = max(d['xmax'] for d in minmax_data)
            ymin = min(d['ymin'] for d in minmax_data)
            ymax = max(d['ymax'] for d in minmax_data)
        else:
            xmin = xmax = ymin = ymax = 0

        if xmin == xmax:
            xmin -= 1
            xmax += 1
        if ymin == ymax:
            ymin -= 1
            ymax += 1

        return {'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax}

    def save(self, filename=None,
             xmin=None, xmax=None, ymin=None, ymax=None,
             figsize=DEFAULT_FIGSIZE, figure=None, sub=None, savenow=True,
             dpi=DEFAULT_DPI, axes=None, axes_labels=None, fontsize=None,
             frame=False, verify=True, aspect_ratio = None,
             gridlines=None, gridlinesstyle=None,
             vgridlinesstyle=None, hgridlinesstyle=None):
        r"""
        Save the graphics to an image file of type: PNG, PS, EPS, SVG, SOBJ,
        depending on the file extension you give the filename.
            Extension types can be: \file{.png}, \file{.ps}, \file{.eps}, \file{.svg},
            and \file{.sobj} (for a \sage object you can load later).

        EXAMPLES:
            sage: c = circle((1,1),1,rgbcolor=(1,0,0))
            sage: c.show(xmin=-1,xmax=3,ymin=-1,ymax=3)

            To correct the apect ratio of certain graphics, it is necessary
            to show with a '\code{figsize}' of square dimensions.

            sage: c.show(figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)

            sage: point((-1,1),pointsize=30, rgbcolor=(1,0,0))

        """
        d = self.get_minmax_data()
        self.axes_range(xmin, xmax, ymin, ymax)
        d.update(self.__axes_range)
        xmin = d['xmin']
        xmax = d['xmax']
        ymin = d['ymin']
        ymax = d['ymax']

        # adjust the figsize in case the user also specifies an aspect ratio
        if aspect_ratio is None:
            aspect_ratio = self.aspect_ratio()
        figsize = adjust_figsize_for_aspect_ratio(figsize, aspect_ratio, xmin=xmin,
                                                  xmax=xmax, ymin=ymin, ymax=ymax)

        global do_verify
        do_verify = verify

        if axes is None:
            axes = self.__show_axes

        from matplotlib.figure import Figure
        if filename is None:
            filename = sage.misc.misc.graphics_filename()
        try:
            ext = os.path.splitext(filename)[1].lower()
        except IndexError:
            raise ValueError, "file extension must be either 'png', 'eps', 'svg' or 'sobj'"

        if ext == '' or ext == '.sobj':
            SageObject.save(self, filename)
            return

        self.fontsize(fontsize)
        self.axes_labels(l=axes_labels)

        if figure is None:
            figure = Figure(figsize)

        #The line below takes away the excessive whitespace around
        #images.  ('figsize' and  'dpi' still work as expected):
        figure.subplots_adjust(left=0.04, bottom=0.04, right=0.96, top=0.96)

        #the incoming subplot instance
        subplot = sub
        if not subplot:
            subplot = figure.add_subplot(111)

        #take away the matplotlib axes:
        subplot.xaxis.set_visible(False)
        subplot.yaxis.set_visible(False)
        subplot.set_frame_on(False)

        #add all the primitives to the subplot
        #check if there are any ContourPlot instances
        #in self._objects, and if so change the axes
        #to be frame axes instead of centered axes
        contour = False
        plotfield = False
        matrixplot = False
        from contour_plot import ContourPlot
        from matrix_plot import MatrixPlot
        from plot_field import PlotField
        for g in self.__objects:
            if isinstance(g, ContourPlot):
                contour = True
            if isinstance(g, PlotField):
                plotfield = True
            if isinstance(g, MatrixPlot):
                matrixplot = True
            g._render_on_subplot(subplot)

        #adjust the xy limits and draw the axes:
        if axes is None:
            axes = self.__show_axes

        #construct an Axes instance, see 'axes.py' for relevant code
        sage_axes = Axes(color=self.__axes_color, fontsize=self.__fontsize,
                         axes_labels=self.__axes_labels,
                         axes_label_color=self.__axes_label_color,
                         tick_label_color=self.__tick_label_color, linewidth=self.__axes_width)

        # construct a GridLines instance, see 'axes.py' for relevant code
        sage_gridlines = GridLines(gridlines=gridlines, gridlinesstyle=gridlinesstyle,
                vgridlinesstyle=vgridlinesstyle, hgridlinesstyle=hgridlinesstyle)

        #adjust the xy limits and draw the axes:
        if not (contour or plotfield or matrixplot): #the plot is a 'regular' plot
            xmin -= 0.1*(xmax-xmin)
            xmax += 0.1*(xmax-xmin)
            ymin -= 0.1*(ymax-ymin)
            ymax += 0.1*(ymax-ymin)
            if frame: #add the frame axes
                axmin, axmax = xmin - 0.04*abs(xmax - xmin), xmax + 0.04*abs(xmax - xmin)
                aymin, aymax = ymin - 0.04*abs(ymax - ymin), ymax + 0.04*abs(ymax - ymin)
                subplot.set_xlim([axmin, axmax])
                subplot.set_ylim([aymin, aymax])
                # draw the grid
                sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, True)
                #add a frame to the plot and possibly 'axes_with_no_ticks'
                sage_axes.add_xy_frame_axes(subplot, xmin, xmax, ymin, ymax,
                                        axes_with_no_ticks=axes)

            elif not frame and axes: #regular plot with regular axes
                # draw the grid
                sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, False)
                # draw the axes
                xmin, xmax, ymin, ymax = sage_axes.add_xy_axes(subplot, xmin, xmax, ymin, ymax)
                subplot.set_xlim(xmin, xmax)
                subplot.set_ylim(ymin, ymax)

            else: #regular plot with no axes
                subplot.set_xlim(xmin, xmax)
                subplot.set_ylim(ymin, ymax)
                # draw the grid
                sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, False)

        elif (contour or plotfield): #contour or field plot in self.__objects, so adjust axes accordingly
            subplot.set_xlim([xmin - 0.05*abs(xmax - xmin), xmax + 0.05*abs(xmax - xmin)])
            subplot.set_ylim([ymin - 0.05*abs(ymax - ymin), ymax + 0.05*abs(ymax - ymin)])
            # draw the grid
            sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, True)
            # draw the axes
            if axes: #axes=True unless user specifies axes=False
                sage_axes.add_xy_frame_axes(subplot, xmin, xmax, ymin, ymax)

        else: #we have a 'matrix_plot' in self.__objects, so adjust axes accordingly
            subplot.set_xlim([xmin - 0.05*abs(xmax - xmin), xmax + 0.05*abs(xmax - xmin)])
            subplot.set_ylim([ymin - 0.05*abs(ymax - ymin), ymax + 0.05*abs(ymax - ymin)])
            # draw the grid
            if gridlines in ["major", "automatic", True]:
                gridlines = [sage.misc.misc.srange(-0.5,xmax+1,1),
                        sage.misc.misc.srange(-0.5,ymax+1,1)]
            sage_gridlines = GridLines(gridlines=gridlines,
                    gridlinesstyle=gridlinesstyle,
                    vgridlinesstyle=vgridlinesstyle,
                    hgridlinesstyle=hgridlinesstyle)
            sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, False)
            # draw the axes
            if axes: #axes=True unless user specifies axes=False
                sage_axes.add_xy_matrix_frame_axes(subplot, xmin, xmax, ymin, ymax)

        # You can output in PNG, PS, EPS, PDF, or SVG format, depending on the file extension.
        # matplotlib looks at the file extension to see what the renderer should be.
        # The default is FigureCanvasAgg for png's because this is by far the most
        # common type of files rendered, like in the Notebook for example.
        # if the file extension is not '.png', then matplotlib will handle it.
        if savenow:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            canvas = FigureCanvasAgg(figure)
            if ext in ['.eps', '.ps', '.pdf']:
                if dpi is None:
                    dpi = 72
            elif ext == '.svg':
                if dpi is None:
                    dpi = 80
            elif ext == '.png':
                if dpi is None:
                    dpi = 100
            else:
                raise ValueError, "file extension must be either 'png', 'ps, 'eps', 'pdf, 'svg' or 'sobj'"
            canvas.print_figure(filename, dpi=dpi)

def xydata_from_point_list(points):
    r"""
    Returns two lists (xdata, ydata), each coerced to a list of
    floats, which correspond to the x-coordinates and the
    y-coordinates of the points.

    The points parameter can be a list of 2-tuples or some object that
    yields a list of one or two numbers.

    This function can potentially be very slow for large point sets.

    """
    if not isinstance(points, (list,tuple)):
        try:
            points = [[float(z) for z in points]]
        except TypeError:
            pass
    elif len(points)==2 and not isinstance(points[0], (list,tuple)):
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
@options(alpha=1, thickness=1, rgbcolor=(0,0,1), plot_points=200,
         adaptive_tolerance=0.01, adaptive_recursion=5, __original_opts=True)
def plot(funcs, *args, **kwds):
    r"""
    Use plot by writing

        \code{plot(X, ...)}

    where $X$ is a \sage object (or list of \sage objects) that either is
    callable and returns numbers that can be coerced to floats, or has
    a plot method that returns a \class{GraphicPrimitive} object.

    Type \code{plot.options} for a dictionary of the default
    options for plots.  You can change this to change
    the defaults for all future plots.  Use \code{plot.reset()}
    to reset to the default options.

    PLOT OPTIONS:
    The plot options are
        plot_points -- (default: 200) the number of points to initially plot
                       before the adaptive refinement.
        adaptive_recursion -- (default: 5) how many levels of recursion to go
                              before giving up when doing adaptive refinement.
                              Setting this to 0 disables adaptive refinement.
        adaptive_tolerance -- (default: 0.01) how large a difference should be
                              before the adaptive refinement code considers it
                              significant.  Use a smaller value for smoother
                              plots and a larger value for coarser plots.  In
                              general, adjust adaptive_recursion before
                              adaptive_tolerance, and consult the code of
                              adaptive_refinement for the details.

        xmin -- starting x value
        xmax -- ending x value
        color -- an rgb-tuple (r,g,b) with each of r,g,b between 0 and 1, or
                 a color name as a string (e.g., 'purple'), or an HTML
                 color such as '\#aaff0b'.

    APPEARANCE OPTIONS:
    The following options affect the appearance of the line through the points
    on the graph of $X$ (these are the same as for the line function):

    INPUT:
        alpha -- How transparent the line is
        thickness -- How thick the line is
        rgbcolor -- The color as an rgb tuple
        hue -- The color given as a hue
        Any MATPLOTLIB line option may also be passed in.  E.g.,
        linestyle -- The style of the line, which is one of
                  '--' (dashed), '-.' (dash dot), '-' (solid),
                  'steps', ':' (dotted)
        marker  -- "'0' (tickleft), '1' (tickright), '2' (tickup), '3' (tickdown),
                   '' (nothing), ' ' (nothing), '+' (plus), ',' (pixel), '.' (point),
                   '1' (tri_down), '3' (tri_left), '2' (tri_up), '4' (tri_right),
                   '<' (triangle_left), '>' (triangle_right), 'None' (nothing),
                   'D' (diamond), 'H' (hexagon2), '_' (hline), '\^' (triangle_up),
                   'd' (thin_diamond), 'h' (hexagon1), 'o' (circle), 'p' (pentagon),
                   's' (square), 'v' (triangle_down), 'x' (x), '|' (vline)"
       markersize -- the size of the marker in points
       markeredgecolor -- the markerfacecolor can be any color arg
       markeredgewidth -- the size of the marker edge in points

    Note that this function does NOT simply sample equally spaced
    points between xmin and xmax.  Instead it computes equally spaced
    points and add small perturbations to them.  This reduces the
    possibility of, e.g., sampling sin only at multiples of $2\pi$,
    which would yield a very misleading graph.

    EXAMPLES:
    We plot the sin function:
        sage: P = plot(sin, (0,10)); print P
        Graphics object consisting of 1 graphics primitive
        sage: len(P)     # number of graphics primitives
        1
        sage: len(P[0])  # how many points were computed (random)
        225
        sage: P          # render

        sage: P = plot(sin, (0,10), plot_points=10); print P
        Graphics object consisting of 1 graphics primitive
        sage: len(P[0])  # random output
        32
        sage: P          # render

    We plot with randomize=False, which makes the initial sample
    points evenly spaced (hence always the same).  Adaptive plotting
    might insert other points, however, unless adaptive_recursion=0.
        sage: p=plot(1, (x,0,3), plot_points=4, randomize=False, adaptive_recursion=0)
        sage: list(p[0])
        [(0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]

    Some colored functions:

        sage: plot(sin, 0, 10, rgbcolor='#ff00ff')
        sage: plot(sin, 0, 10, rgbcolor='purple')

    We plot several functions together by passing a list
    of functions as input:
        sage: plot([sin(n*x) for n in [1..4]], (0, pi))


    The function $\sin(1/x)$ wiggles wildly near $0$.  Sage adapts
    to this and plots extra points near the origin.
        sage: plot(sin(1/x), (x, -1, 1))

    Note that the independent variable may be omitted if there is no
    ambiguity:
        sage: plot(sin(1/x), (-1, 1))

    The algorithm used to insert extra points is actually pretty simple. On
    the picture drawn by the lines below:
        sage: p = plot(x^2, (-0.5, 1.4)) + line([(0,0), (1,1)], rgbcolor='green')
        sage: p += line([(0.5, 0.5), (0.5, 0.5^2)], rgbcolor='purple')
        sage: p += point(((0, 0), (0.5, 0.5), (0.5, 0.5^2), (1, 1)), rgbcolor='red', pointsize=20)
        sage: p += text('A', (-0.05, 0.1), rgbcolor='red')
        sage: p += text('B', (1.01, 1.1), rgbcolor='red')
        sage: p += text('C', (0.48, 0.57), rgbcolor='red')
        sage: p += text('D', (0.53, 0.18), rgbcolor='red')
        sage: p.show(axes=False, xmin=-0.5, xmax=1.4, ymin=0, ymax=2)

    You have the function (in blue) and its approximation (in green) passing
    through the points A and B. The algorithm finds the midpoint C of AB and
    computes the distance between C and D. The point D is added to the curve if
    it exceeds the (nonzero) adaptive_tolerance threshold. If D is added to
    the curve, then the algorithm is applied recursively to the points A and D,
    and D and B. It is repeated adaptive_recursion times (10, by default).

    The actual sample points are slightly randomized, so the above
    plots may look slightly different each time you draw them.

    We draw the graph of an elliptic curve as the union
    of graphs of 2 functions.
        sage: def h1(x): return abs(sqrt(x^3  - 1))
        sage: def h2(x): return -abs(sqrt(x^3  - 1))
        sage: P = plot([h1, h2], 1,4)
        sage: P          # show the result

    We can also directly plot the elliptic curve:
        sage: E = EllipticCurve([0,-1])
        sage: plot(E, (1, 4), rgbcolor=hue(0.6))

    We can change the line style to one of '--' (dashed), '-.' (dash dot),
    '-' (solid), 'steps', ':' (dotted):
        sage: plot(sin(x), 0, 10, linestyle='-.')

    Sage currently ignores points that cannot be evaluated
        sage: plot(-x*log(x), (x,0,1))  # this works fine since the failed endpoint is just skipped.

    This prints out a warning and plots where it can (we turn off the warning by setting
    the verbose mode temporarily to -1.)
        sage: set_verbose(-1)
        sage: plot(x^(1/3), (x,-1,1))
        sage: set_verbose(0)

    To plot the negative real cube root, use something like the following.
        sage: plot(lambda x : RR(x).nth_root(3), (x,-1, 1))

    TESTS:
    We do not randomize the endpoints:
        sage: p = plot(x, (x,-1,1))
        sage: p[0].xdata[0] == -1
        True
        sage: p[0].xdata[-1] == 1
        True

    We check to make sure that the x/y min/max data get set correctly
    when there are multiple functions.

        sage: d = plot([sin(x), cos(x)], 100, 120).get_minmax_data()
        sage: d['xmin']
        100.0
        sage: d['xmax']
        120.0

    We check to handle cases where the function gets evaluated at a point
    which causes an 'inf' or '-inf' result to be produced.
        sage: p = plot(1/x, 0, 1)
        sage: p = plot(-1/x, 0, 1)
    """
    original_opts = kwds.pop('__original_opts', {})
    do_show = kwds.pop('show',False)
    if hasattr(funcs, 'plot'):
        G = funcs.plot(*args, **original_opts)
    # if we are using the generic plotting method
    else:
        n = len(args)
        # if there are no extra args, pick some silly default
        if n == 0:
            G = _plot(funcs, (-1, 1), *args, **kwds)
        # if there is one extra arg, then it had better be a tuple
        elif n == 1:
            G = _plot(funcs, *args, **kwds)
        elif n == 2:
        # if there are two extra args, then pull them out and pass them as a tuple
            xmin = args[0]
            xmax = args[1]
            args = args[2:]
            G = _plot(funcs, (xmin, xmax), *args, **kwds)
        else:
            sage.misc.misc.verbose("there were %s extra arguments (besides %s)" % (n, funcs), level=0)
    if kwds.has_key('xmin') and kwds.has_key('xmax'):
        xmin = kwds['xmin']
        xmax = kwds['xmax']
        del kwds['xmin']
        del kwds['xmax']
        G = _plot(funcs, (xmin, xmax), *args, **kwds)
    if do_show:
        G.show()
    return G

def _plot(funcs, xrange, parametric=False,
              polar=False, label='', randomize=True, **options):
    if not is_fast_float(funcs):
        funcs =  fast_float(funcs)

    #parametric_plot will be a list or tuple of two functions (f,g)
    #and will plotted as (f(x), g(x)) for all x in the given range
    if parametric:
        f, g = funcs
    #or we have only a single function to be plotted:
    else:
        f = funcs

    plot_points = int(options.pop('plot_points'))
    x, data = var_and_list_of_values(xrange, plot_points)
    xmin = data[0]
    xmax = data[-1]

    #check to see if funcs is a list of functions that will
    #be all plotted together.
    if isinstance(funcs, (list, tuple)) and not parametric:
        return reduce(operator.add, (plot(f, (xmin, xmax), polar=polar, **options) for f in funcs))

    delta = float(xmax-xmin) / plot_points

    random = current_randstate().python_random().random
    exceptions = 0; msg=''
    exception_indices = []
    for i in range(len(data)):
        xi = data[i]
        # Slightly randomize the interior sample points if
        # randomize is true
        if randomize and i > 0 and i < plot_points-1:
            xi += delta*(random() - 0.5)

        try:
            data[i] = (float(xi), float(f(xi)))
            if str(data[i][1]) in ['nan', 'NaN', 'inf', '-inf']:
                sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, x),1)
                exceptions += 1
                exception_indices.append(i)
        except (ZeroDivisionError, TypeError, ValueError, OverflowError), msg:
            sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, x),1)

            if i == 0:
                for j in range(1, 99):
                    xj = xi + delta*j/100.0
                    try:
                        data[i] = (float(xj), float(f(xj)))
                        # nan != nan
                        if data[i][1] != data[i][1]:
                            continue
                        break
                    except (ZeroDivisionError, TypeError, ValueError, OverflowError), msg:
                        pass
                else:
                    exceptions += 1
                    exception_indices.append(i)
            elif i == plot_points-1:
                for j in range(1, 99):
                    xj = xi - delta*j/100.0
                    try:
                        data[i] = (float(xj), float(f(xj)))
                        # nan != nan
                        if data[i][1] != data[i][1]:
                            continue
                        break
                    except (ZeroDivisionError, TypeError, ValueError, OverflowError), msg:
                        pass
                else:
                    exceptions += 1
                    exception_indices.append(i)
            else:
                exceptions += 1
                exception_indices.append(i)

            exceptions += 1
            exception_indices.append(i)


    data = [data[i] for i in range(len(data)) if i not in exception_indices]

    # adaptive refinement
    i, j = 0, 0
    adaptive_tolerance = delta * float(options.pop('adaptive_tolerance'))
    adaptive_recursion = int(options.pop('adaptive_recursion'))

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
    if parametric:
        data = [(fdata, g(x)) for x, fdata in data]
    if polar:
        data = [(y*cos(x), y*sin(x)) for x, y in data]

    from sage.plot.all import line, text
    G = line(data, **options)

    # Label?
    if label:
        label = '  '+str(label)
        G += text(label, data[-1], horizontal_alignment='left',
                  vertical_alignment='center')

    return G




########## misc functions ###################

def parametric_plot(funcs, tmin, tmax, **kwargs):
    r"""
    \code{parametric_plot} takes two or three functions as a list or a
    tuple and makes a plot with the first function giving the $x$
    coordinates, the second function giving the $y$ coordinates, and the
    third function (if present) giving the $z$ coordinates.

    INPUT:
        funcs -- 2 or 3-tuple of functions
        tmin -- start value of t
        tmax -- end value of t
        other options -- passed to plot.

    EXAMPLES:
    We draw some 2d parametric plots:
        sage: t = var('t')
        sage: parametric_plot( (sin(t), sin(2*t)), 0, 2*pi, rgbcolor=hue(0.6) )
        sage: parametric_plot((1, t), 0, 4)
        sage: parametric_plot((t, t^2), -4, 4)

    We draw a 3d parametric plot:
        sage: parametric_plot3d( (5*cos(x), 5*sin(x), x), (-12, 12), plot_points=150, color="red")

    TESTS:
        sage: parametric_plot((x, t^2), -4, 4)
        Traceback (most recent call last):
        ...
        ValueError: there cannot be more than one free variable in funcs

        sage: parametric_plot((1, x+t), -4, 4)
        Traceback (most recent call last):
        ...
        ValueError: there cannot be more than one free variable in funcs

    """
    if len(funcs) == 3:
        raise ValueError, "use parametric_plot3d for parametric plots in 3d dimensions."
    elif len(funcs) != 2:
        raise ValueError, "parametric plots only implemented in 2 and 3 dimensions."
    else:
        vars = []
        f,g = funcs
        if hasattr(f, 'variables'):
            vars += list(f.variables())
        if hasattr(g, 'variables'):
            vars += list(g.variables())
        vars = [str(v) for v in vars]

        from sage.misc.misc import uniq
        if len(uniq(vars)) > 1:
            raise ValueError, "there cannot be more than one free variable in funcs"

    return plot(funcs, tmin, tmax, parametric=True, **kwargs)

def polar_plot(funcs, xmin, xmax, **kwargs):
    r"""
    \code{polar_plot} takes a single function or a list or tuple of functions
    and plots them parametrically in the given range.

    EXAMPLES:
    Here is a blue 8-leaved petal:
        sage: polar_plot(sin(5*x)^2, 0, 2*pi, rgbcolor=hue(0.6))

    A red figure-8:
        sage: polar_plot(abs(sqrt(1 - sin(x)^2)), 0, 2*pi, rgbcolor=hue(1.0))

    A green limacon of Pascal:
        sage: polar_plot(2 + 2*cos(x), 0, 2*pi, rgbcolor=hue(0.3))

    """
    return plot(funcs, xmin, xmax, polar=True, **kwargs)

def list_plot(data, plotjoined=False, **kwargs):
    r"""
    \code{list_plot} takes a single list of data, in which case it forms a
    list of tuples $(i,di)$ where $i$ goes from 0 to ${\rm len}(data)-1$ and $di$ is
    the $i$th data value, and puts points at those tuple values.

    \code{list_plot} also takes a list of tuples $(dxi, dyi)$ where $dxi$ is the
    $i$th data representing the $x$-value, and $dyi$ is the $i$th $y$-value.  If
    \code{plotjoined=True}, then a line spanning all the data is drawn
    instead.

    EXAMPLES:
        sage: list_plot([i^2 for i in range(5)])

    Here are a bunch of random red points:
        sage: r = [(random(),random()) for _ in range(20)]
        sage: list_plot(r,rgbcolor=(1,0,0))

    This gives all the random points joined in a purple line:
        sage: list_plot(r, plotjoined=True, rgbcolor=(1,0,1))

    TESTS:
    We check to see that the x/y min/max data are set correctly.
        sage: d = list_plot([(100,100), (120, 120)]).get_minmax_data()
        sage: d['xmin']
        100.0
        sage: d['ymin']
        100.0

    """
    from sage.plot.all import line, point
    if not isinstance(data[0], (list, tuple)):
        data = zip(range(len(data)),data)
    if plotjoined:
        P = line(data, **kwargs)
    else:
        P = point(data, **kwargs)
    return P


def to_float_list(v):
    """
    Given a list or tuple or iterable v, coerce each element of v to a
    float and make a list out of the result.

    EXAMPLES:
        sage: from sage.plot.plot import to_float_list
        sage: to_float_list([1,1/2,3])
        [1.0, 0.5, 3.0]
    """
    return [float(x) for x in v]


def hue(h, s=1, v=1):
    """
      hue(h,s=1,v=1) where 'h' stands for hue,
      's' stands for saturation, 'v' stands for value.
      hue returns a tuple of rgb intensities (r, g, b)
      All values are in the range 0 to 1.

      INPUT:
         h, s, v -- real numbers between 0 and 1.  Note that
                    if any are not in this range they are automatically
                    normalized to be in this range by reducing them
                    modulo 1.
      OUTPUT:
         A valid RGB tuple.

      EXAMPLES:
        sage: hue(0.6)
        (0.0, 0.40000000000000036, 1.0)

      hue is an easy way of getting a broader
      range of colors for graphics

        sage: plot(sin, -2, 2, rgbcolor=hue(0.6))

    """
    h = float(h); s = float(s); v = float(v)
    if h != 1:
        h = modf(h)[0]
        if h < 0:
            h += 1
    if s != 1:
        s = modf(s)[0]
        if s < 0:
            s += 1
    if v != 1:
        v = modf(v)[0]
        if v < 0:
            v += 1
    c = hsv_to_rgb(h, s, v)
    return (float(c[0]), float(c[1]), float(c[2]))

class GraphicsArray(SageObject):
    """
    GraphicsArray takes a ($m$ x $n$) list of lists of graphics
    objects and plots them all on one canvas.
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
        self._figsize = DEFAULT_FIGSIZE

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

    def _render(self, filename, dpi=None, figsize=DEFAULT_FIGSIZE, axes=None, **args):
        r"""
        \code{render} loops over all graphics objects
        in the array and adds them to the subplot.
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
            g.save(filename, dpi=dpi, figure=figure, sub=subplot,
                   savenow = (i==dims), verify=do_verify,
                   axes = axes,
                   **args)#only save if i==dims.

    def save(self, filename=None, dpi=DEFAULT_DPI, figsize=DEFAULT_FIGSIZE,
             axes = None, **args):
        """
        save the \code{graphics_array} to
            (for now) a png called 'filename'.
        """
        if (figsize != DEFAULT_FIGSIZE): self.__set_figsize__(figsize)
        self._render(filename, dpi=dpi, figsize=self._figsize, axes = axes, **args)

    def show(self, filename=None, dpi=DEFAULT_DPI, figsize=DEFAULT_FIGSIZE,
             axes = None, **args):
        r"""
        Show this graphics array using the default viewer.

        OPTIONAL INPUT:
            filename -- (default: None) string
            dpi -- dots per inch
            figsize -- [width, height] (same for square aspect)
            axes -- (default: True)
            fontsize -- positive integer
            frame -- (default: False) draw a frame around the image

        EXAMPLES:
        This draws a graphics array with four trig plots and no axes
        in any of the plots.
            sage: G = graphics_array([[plot(sin), plot(cos)], [plot(tan), plot(sec)]])
            sage: G.show(axes=False)
        """
        if (figsize != DEFAULT_FIGSIZE): self.__set_figsize__(figsize)
        if DOCTEST_MODE:
            self.save(sage.misc.misc.SAGE_TMP + '/test.png',
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
    \code{graphics_array} take a list of lists (or tuples)
    of graphics objects and plots them all on one canvas (single plot).

    INPUT:
         array -- a list of lists or tuples
         n, m -- (optional) integers -- if n and m are given then
                 the input array is flattened and turned into an
                 n x m array, with blank graphics objects padded
                 at the end, if necessary.

    EXAMPLE:
    Make some plots of $\sin$ functions:

        sage: f(x) = sin(x)
        sage: g(x) = sin(2*x)
        sage: h(x) = sin(4*x)
        sage: p1 = plot(f,-2*pi,2*pi,rgbcolor=hue(0.5))
        sage: p2 = plot(g,-2*pi,2*pi,rgbcolor=hue(0.9))
        sage: p3 = parametric_plot((f,g),0,2*pi,rgbcolor=hue(0.6))
        sage: p4 = parametric_plot((f,h),0,2*pi,rgbcolor=hue(1.0))

    Now make a graphics array out of the plots;
    then you can type either: \code{ga.show()} or \code{ga.save()}.

        sage: graphics_array(((p1,p2),(p3,p4)))

    Here we give only one row:
        sage: p1 = plot(sin,-4,4)
        sage: p2 = plot(cos,-4,4)
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

def float_to_html(r,g,b):
    """
    This is a function to present tuples of RGB floats as HTML-happy hex
    for matplotlib. This may not seem necessary, but there are some odd
    cases where matplotlib is just plain schizophrenic--for an example, do

    EXAMPLES:
        sage: vertex_colors = {(1.0, 0.8571428571428571, 0.0): [4, 5, 6], (0.28571428571428559, 0.0, 1.0): [14, 15, 16], (1.0, 0.0, 0.0): [0, 1, 2, 3], (0.0, 0.57142857142857162, 1.0): [12, 13], (1.0, 0.0, 0.85714285714285676): [17, 18, 19], (0.0, 1.0, 0.57142857142857162): [10, 11], (0.28571428571428581, 1.0, 0.0): [7, 8, 9]}
        sage: graphs.DodecahedralGraph().show(vertex_colors=vertex_colors)

    Notice how the colors don't respect the partition at all.....
    """ # TODO: figure out WTF
    from sage.rings.integer import Integer
    from math import floor
    rr = Integer(int(floor(r*255))).str(base=16)
    gg = Integer(int(floor(g*255))).str(base=16)
    bb = Integer(int(floor(b*255))).str(base=16)
    rr = '0'*(2-len(rr)) + rr
    gg = '0'*(2-len(gg)) + gg
    bb = '0'*(2-len(bb)) + bb
    return '#' + rr\
               + gg\
               + bb

def rainbow(n, format='hex'):
    """
    Given an integer $n$, returns a list of colors, represented in HTML hex,
    that changes smoothly in hue from one end of the spectrum to the other.
    Written in order to easily represent vertex partitions on graphs.

    AUTHOR: Robert L. Miller

    EXAMPLE:
        sage: from sage.plot.plot import rainbow
        sage: rainbow(7)
        ['#ff0000', '#ffda00', '#48ff00', '#00ff91', '#0091ff', '#4800ff', '#ff00da']
        sage: rainbow(7, 'rgbtuple')
        [(1.0, 0.0, 0.0), (1.0, 0.8571428571428571, 0.0), (0.28571428571428581, 1.0, 0.0), (0.0, 1.0, 0.57142857142857162), (0.0, 0.57142857142857162, 1.0), (0.28571428571428559, 0.0, 1.0), (1.0, 0.0, 0.85714285714285676)]

    """
    from math import floor
    R = []
    if format == 'hex':
        for i in range(n):
            r = 6*float(i)/n
            h = floor(r)
            r = float(r - h)
            if h == 0:#RED
                R.append(float_to_html(1.,r,0.))
            elif h == 1:
                R.append(float_to_html(1. - r,1.,0.))
            elif h == 2:#GREEN
                R.append(float_to_html(0.,1.,r))
            elif h == 3:
                R.append(float_to_html(0.,1. - r,1.))
            elif h == 4:#BLUE
                R.append(float_to_html(r,0.,1.))
            elif h == 5:
                R.append(float_to_html(1.,0.,1. - r))
    elif format == 'rgbtuple':
        for i in range(n):
            r = 6*float(i)/n
            h = floor(r)
            r = float(r - h)
            if h == 0:#RED
                R.append((1.,r,0.))
            elif h == 1:
                R.append((1. - r,1.,0.))
            elif h == 2:#GREEN
                R.append((0.,1.,r))
            elif h == 3:
                R.append((0.,1. - r,1.))
            elif h == 4:#BLUE
                R.append((r,0.,1.))
            elif h == 5:
                R.append((1.,0.,1. - r))
    return R

def var_and_list_of_values(v, plot_points):
    """
    INPUT:
        v -- (v0, v1) or (var, v0, v1); if the former return
             the range of values between v0 and v1 taking
             plot_points steps; if var is given, also return var.
        plot_points -- integer >= 2 (the endpoints)

    OUTPUT:
        var -- a variable or None
        list -- a list of floats

    EXAMPLES:
        sage: from sage.plot.plot import var_and_list_of_values
        sage: var_and_list_of_values((var('theta'), 2, 5),  5)
        (theta, [2.0, 2.75, 3.5, 4.25, 5.0])
        sage: var_and_list_of_values((2, 5),  5)
        (None, [2.0, 2.75, 3.5, 4.25, 5.0])
        sage: var_and_list_of_values((var('theta'), 2, 5),  2)
        (theta, [2.0, 5.0])
        sage: var_and_list_of_values((2, 5),  2)
        (None, [2.0, 5.0])
    """
    plot_points = int(plot_points)
    if plot_points < 2:
        raise ValueError, "plot_points must be positive"
    if not isinstance(v, (tuple, list)):
        raise TypeError, "v must be a tuple or list"
    if len(v) == 3:
        var = v[0]
        a, b = v[1], v[2]
    elif len(v) == 2:
        var = None
        a, b = v
    else:
        raise ValueError, "parametric value range must be a list of 2 or 3-tuple."

    a = float(a)
    b = float(b)
    if plot_points == 2:
        return var, [a, b]
    else:
        step = (b-a)/float(plot_points-1)
        values = [a + step*i for i in xrange(plot_points)]
        return var, values


def adjust_figsize_for_aspect_ratio(figsize, aspect_ratio, xmin, xmax, ymin, ymax):
    """
    Adjust the figsize in case the user also specifies an aspect ratio.

    INPUTS:
        figsize -- a sequence of two positive real numbers
        aspect_ratio -- a positive real number
        xmin, xmax, ymin, ymax -- real numbers

    EXAMPLES:
    This function is used mainly internally by plotting code so we explicitly import it:
        sage: from sage.plot.plot import adjust_figsize_for_aspect_ratio

    This returns (5,5), since the requested aspect ratio is 1 and the
    x and y ranges are the same, so that's the right size rendered
    image to produce a 1:1 ratio internally.  5 is used instead of 3
    since the image size is always adjusted to the larger of the
    figsize dimensions.

        sage: adjust_figsize_for_aspect_ratio([3,5], 1, 0, 2, 0, 2)
        (5, 5)

    Here we give a scalar figsize, which is automatically converted to
    the figsize \code{(figsize, figsize/golden_ratio)}.
        sage: adjust_figsize_for_aspect_ratio(3, 1, 0, 2, 0, 2)
        (3, 3)

    Here we omit the aspect ratio so the figsize is just returned.
        sage: adjust_figsize_for_aspect_ratio([5,6], None, 0, 2, 0, 2)
        [5, 6]

    Here we have an aspect ratio of 2, and since the x and y ranges are
    the same the returned figsize is twice as wide as tall:
        sage: adjust_figsize_for_aspect_ratio([3,5], 2, 0, 2, 0, 2)
        (5, 5/2)

    Here the x range is rather large, so to get an aspect ratio where circles
    look twice as wide as they are tall, we have to shrink the y size
    of the image.
        sage: adjust_figsize_for_aspect_ratio([3,5], 2, 0, 10, 0, 2)
        (5, 1/2)
    """
    if not isinstance(figsize, (list, tuple)):
        figsize = [figsize, figsize * 0.618033988749895]   # 1/golden_ratio
    if aspect_ratio is None:
        return figsize
    # We find a number r such that (ymax-ymin)*r / (xmax-xmin) = aspect_ratio:
    r = max(aspect_ratio * (xmax - xmin)/(ymax-ymin), 0.001)
    mx = max(figsize)
    f = (figsize[0]*r, figsize[0])
    s = min((mx/f[0], mx/f[1]))
    return f[0]*s, f[1]*s



def setup_for_eval_on_grid(v, xrange, yrange, plot_points):
    """
    INPUT:
        v -- a list of functions
        xrange -- 2 or 3 tuple (if 3, first is a variable)
        yrange -- 2 or 3 tuple
        plot_points -- a positive integer

    OUTPUT:
        g -- tuple of fast callable functions
        xstep -- step size in xdirection
        ystep -- step size in ydirection
        xrange -- tuple of 2 floats
        yrange -- tuple of 2 floats

    EXAMPLES:
        sage: x,y = var('x,y')
        sage: sage.plot.plot.setup_for_eval_on_grid([x^2 + y^2], (x,0,5), (y,0,pi), 10)
        ([<sage.ext.fast_eval.FastDoubleFunc object at ...>],
         0.5,
         0.31415926535897931,
         (0.0, 5.0),
         (0.0, 3.1415926535897931))
    """
    if len(xrange) == 3:
        xvar = xrange[0]
        xrange = xrange[1:]
        yvar = yrange[0]
        yrange = yrange[1:]
    else:
        xvar = None
    xrange = tuple([float(z) for z in xrange])
    yrange = tuple([float(z) for z in yrange])
    plot_points = int(plot_points)
    if plot_points <= 0:
        plot_points = 1
    xstep = abs(xrange[0] - xrange[1])/plot_points
    ystep = abs(yrange[0] - yrange[1])/plot_points

    g = []
    for f in v:
        if isinstance(f, types.FunctionType):
            g.append(f)
        else:
            # This code can be refactored at some point out of plot3d.
            from sage.plot.plot3d.parametric_plot3d import adapt_to_callable
            if xvar is None:
                k, _ = adapt_to_callable([f], 2)
                g.append(k[0])
            else:
                g.append(fast_float(f, str(xvar), str(yvar)))

    return g, xstep, ystep, xrange, yrange


def minmax_data(xdata, ydata, dict=False):
    """
    Returns the minimums and maximums of xdata and ydata.

    If dict is False, then minmax_data returns the tuple
    (xmin, xmax, ymin, ymax); otherwise, it returns a dictionary
    whose keys are 'xmin', 'xmax', 'ymin', and 'ymax' and whose
    values are the corresponding values.

    EXAMPLES:
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
    The adaptive refinement algorithm for plotting a function f. See the
    docstring for plot for a description of the algorithm.

    INPUT:
        f -- a function of one variable
        p1, p2 -- two points to refine between
        adaptive_recursion -- (default: 5) how many levels of recursion to go
                              before giving up when doing adaptive refinement.
                              Setting this to 0 disables adaptive refinement.
        adaptive_tolerance -- (default: 0.01) how large a difference should be
                              before the adaptive refinement code considers
                              it significant.  See the documentation for
                              plot() for more information.

    OUTPUT:
        list -- a list of points to insert between p1 and p2 to get
                a better linear approximation between them

    TESTS:
        sage: from sage.plot.plot import adaptive_refinement
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01, adaptive_recursion=0)
        []
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01)
        [(0.125*pi, 0.38268343236508978), (0.1875*pi, 0.55557023301960218), (0.25*pi, 0.70710678118654746), (0.3125*pi, 0.831469612302545...), (0.375*pi, 0.92387953251128674), (0.4375*pi, 0.98078528040323043), (0.5*pi, 1.0), (0.5625*pi, 0.98078528040323043), (0.625*pi, 0.92387953251128674), (0.6875*pi, 0.83146961230254546), (0.75*pi, 0.70710678118654757), (0.8125*pi, 0.55557023301960218), (0.875*pi, 0.3826834323650898...)]

    This shows that lowering adaptive_tolerance and raising
    adaptive_recursion both increase the number of subdivision points:

        sage: x = var('x')
        sage: f = sin(1/x)
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
    try:
        y = float(f(x))
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
