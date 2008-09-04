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

from misc import rgbcolor, Color

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

    def __init__(self, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Create a new empty Graphics objects with all the defaults.

        EXAMPLES:
            sage: G = Graphics()
        """
        self.__xmin = -1 if xmin is None else xmin
        self.__xmax = 1 if xmax is None else xmax
        self.__ymin = -1 if ymin is None else ymin
        self.__ymax = 1 if ymax is None else ymax
        self.__aspect_ratio = None
        self.__fontsize = 10
        self.__show_axes = True
        self.__axes_color = (0, 0, 0)
        self.__axes_label_color = (0, 0, 0)
        self.__tick_label_color = (0, 0, 0)
        self.__axes_width = 0.8
        self.__objects = []

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
            sage: L.xmin(), L.xmax(), L.ymin(), L.ymax()
            (-1.0, 20.0, 0.0, 2.0)
        """
        self.xmin(xmin)
        self.xmax(xmax)
        self.ymin(ymin)
        self.ymax(ymax)

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

    def xmax(self, new=None):
        """
        EXAMPLES:
            sage: G = Graphics(); print G
            Graphics object consisting of 0 graphics primitives
            sage: G.xmax()
            1
            sage: G.xmax(2)
            2.0
            sage: G.xmax()
            2.0
        """
        if new is None:
            return self.__xmax
        new = float(new)
        self.__xmax = new
        return new

    def xmin(self, new=None):
        """
        EXAMPLES:
            sage: G = Graphics(); print G
            Graphics object consisting of 0 graphics primitives
            sage: G.xmin()
            -1
            sage: G.xmax(2)
            2.0
            sage: G.xmax()
            2.0
        """
        if new is None:
            return self.__xmin
        new = float(new)
        self.__xmin = new
        return new

    def ymax(self, new=None):
        """
        EXAMPLES:
            sage: G = Graphics(); print G
            Graphics object consisting of 0 graphics primitives
            sage: G.ymax()
            1
            sage: G.ymax(2)
            2.0
            sage: G.ymax()
            2.0
        """
        if new is None:
            return self.__ymax
        new = float(new)
        self.__ymax = new
        return new

    def ymin(self, new=None):
        """
        EXAMPLES:
            sage: G = Graphics(); print G
            Graphics object consisting of 0 graphics primitives
            sage: G.ymin()
            -1
            sage: G.ymin(2)
            2.0
            sage: G.ymin()
            2.0
        """
        if new is None:
            return self.__ymin
        new = float(new)
        self.__ymin = new
        return new

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
        return self.__objects[int(i)]

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
        g.__xmin = min(self.__xmin, other.__xmin)
        g.__xmax = max(self.__xmax, other.__xmax)
        g.__ymin = min(self.__ymin, other.__ymin)
        g.__ymax = max(self.__ymax, other.__ymax)
        g.__objects = self.__objects + other.__objects
        g.__aspect_ratio = max(self.__aspect_ratio, other.__aspect_ratio)
        return g

    def add_primitive(self, primitive):
        """
        Adds a primitive to this graphics object.
        """
        self.__objects.append(primitive)

    def _arrow(self, xtail, ytail, xhead, yhead, options):
        """
        Add an arrow with given bounding box to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            xtail, ytail, xhead, yhead -- start and stop point of arrow
            options -- dictionary

        EXAMPLES:
        This will display a bold green arrow going up and to the right.
            sage: S = circle((0,0), 2)
            sage: S._arrow(0,0,5,5, {'width':0.2, 'rgbcolor':(0,1,0)}); S
        """
        self.__objects.append(GraphicPrimitive_Arrow(xtail, ytail, xhead, yhead, options))
        self._extend_axes(xtail, xhead, ytail, yhead)

    def _bar_chart(self, ind, datalist, xrange, yrange, options):
        """
        Add a bar chart to this graphics objects.

        (For internal use -- you should just use addition.)

        INPUT:
            ind -- index list
            datalist -- list of values for each element of the index list.
            xrange -- pair (xmin, xmax) of floats
            yrange -- pair (ymin, ymax) of floats
            options -- dictionary of options

        EXAMPLES:
            sage: S = circle((0,0), 2)
            sage: S._bar_chart(range(4), [1,3,2,0], (0,4), (0,3), {'width': 0.5, 'rgbcolor': (0, 0, 1)})
            sage: S
        """
        self.__objects.append(GraphicPrimitive_BarChart(ind, datalist, options))
        self._extend_axes(xrange[0], xrange[1], yrange[0], yrange[1])

    def _circle(self, x, y, r, options):
        """
        Add a circle to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            x -- float; x coordinate of center of circle
            y -- float; y coordinate of center of circle
            r -- float; radius of circle
            options -- dictionary of options

        EXAMPLES:
        We inscribe a circle in another circle:
            sage: S = circle((0,0), 2)
            sage: S._circle(1,0,1, {'rgbcolor':(0.5,0,1), 'thickness':1, 'fill':True, 'alpha':1})
            sage: S.show(aspect_ratio=1, axes=False)
        """
        self.__objects.append(GraphicPrimitive_Circle(x, y, r, options))
        self._extend_axes(x+r, x-r, y+r, y-r)

    def _contour_plot(self, xy_data_array, xrange, yrange, options):
        """
        Add a countor plot to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            xy_data_array --
            xrange, yrange --
            options -- dictionary of options
        """
        self.__xmin = xrange[0]
        self.__xmax = xrange[1]
        self.__ymin = yrange[0]
        self.__ymax = yrange[1]
        self.__objects.append(GraphicPrimitive_ContourPlot(xy_data_array, xrange, yrange, options))

    def _disk(self, point, r, angle, options):
        """
        Add a disk to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            point --
            r --
            angle --
            options -- dictionary of options
        """
        xmin = point[0] - r
        xmax = point[0] + r
        ymin = point[1] - r
        ymax = point[1] + r
        self.__objects.append(GraphicPrimitive_Disk(point, r, angle, options))
        self._extend_axes(xmin, xmax, ymin, ymax)


    def _matrix_plot(self, xy_data_array, xrange, yrange, options):
        """
        Add a matrix plot to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            xy_data_array --
            xrange, yrange --
            options -- dictionary of options
        """
        self.__xmin = xrange[0]
        self.__xmax = xrange[1]
        self.__ymin = yrange[0]
        self.__ymax = yrange[1]
        self.__objects.append(GraphicPrimitive_MatrixPlot(xy_data_array, xrange, yrange, options))
        #self._extend_axes(xrange[0], xrange[1], yrange[0], yrange[1])

    def _plot_field(self, xpos_array, ypos_array, xvec_array, yvec_array, xrange, yrange, options):
        """
        Add a vector field plot to this graphics object.

        INPUT:

            xrange, yrange --
            options -- dictionary of options
        """
        self.__xmin = xrange[0]
        self.__xmax = xrange[1]
        self.__ymin = yrange[0]
        self.__ymax = yrange[1]
        self.__objects.append(GraphicPrimitive_PlotField(xpos_array, ypos_array, xvec_array, yvec_array, options))

    def _point(self, xdata, ydata, options, extend_axes=True):
        """
        Add a plot of a point or list of points to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            xy_data_array --
            xrange, yrange --
            options -- dictionary of options
        """
        self.__objects.append(GraphicPrimitive_Point(xdata, ydata, options))
        if extend_axes:
            self._extend_axes(*minmax_data(xdata, ydata))

    def _polygon(self, xdata, ydata, options, extend_axes=True):
        """
        Add a plot of a polygon to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            xdata -- x coordinates of vertices of the polygon
            ydata -- y coordinates of vertices of the polygon
            options -- dictionary of options
        """
        self.__objects.append(GraphicPrimitive_Polygon(xdata, ydata, options))
        if extend_axes:
            self._extend_axes(*minmax_data(xdata, ydata))

    def _text(self, string, point, options):
        """
        Add a string of text to this graphics object.

        (For internal use -- you should just use addition.)

        INPUT:
            string -- a string
            point -- a 2-tuple (x,y) of floats
            options -- dictionary of options
        """
        self.__objects.append(GraphicPrimitive_Text(string, point, options))
        xpad = 0.2*abs(point[0])
        ypad = 0.2*abs(point[1])
        self._extend_axes(point[0] - xpad, point[0] + xpad, point[1] - ypad, point[1] + ypad)

    def _extend_x_axis(self, x):
        """
        Extend the x axis range so that it contains x.

        EXAMPLES:
            sage: S = circle((0,0), 2)
            sage: S.xmin(), S.xmax()
            (-2.0, 2.0)
            sage: S._extend_x_axis(1)
            sage: S.xmin(), S.xmax()
            (-2.0, 2.0)
            sage: S._extend_x_axis(-5)
            sage: S.xmin(), S.xmax()
            (-5, 2.0)
            sage: S._extend_x_axis(5)
            sage: S.xmin(), S.xmax()
            (-5, 5)
        """
        xmin = self.__xmin
        xmax = self.__xmax
        if xmin is None or x < xmin:
            self.__xmin = x
        elif xmax is None or x > xmax:
            self.__xmax = x

    def _extend_y_axis(self, y):
        """
        Extend the y axis range so that it contains y.

        EXAMPLES:
            sage: S = circle((0,0), 2)
            sage: S.ymin(), S.ymax()
            (-2.0, 2.0)
            sage: S._extend_y_axis(1)
            sage: S.ymin(), S.ymax()
            (-2.0, 2.0)
            sage: S._extend_y_axis(-5)
            sage: S.ymin(), S.ymax()
            (-5, 2.0)
            sage: S._extend_y_axis(5)
            sage: S.ymin(), S.ymax()
            (-5, 5)
        """
        ymin = self.__ymin
        ymax = self.__ymax
        if ymin is None or y < ymin:
            self.__ymin = y
        elif ymax is None or y > ymax:
            self.__ymax = y

    def _extend_axes(self, xmin, xmax, ymin, ymax):
        """
        Extend both the x and y axis so that the x-axis contains both
        xmin and xmax, and the y-axis contains by ymin and ymax.

        EXAMPLES:
            sage: S = circle((0,0), 2)
            sage: S._extend_axes(-2, 3, -5, 1)
            sage: S.xmin(), S.xmax(), S.ymin(), S.ymax()
            (-2.0, 3, -5, 2.0)
        """
        self._extend_x_axis(xmin)
        self._extend_x_axis(xmax)
        self._extend_y_axis(ymin)
        self._extend_y_axis(ymax)

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

    def _prepare_axes(self, xmin, xmax, ymin, ymax):
        """
        Perform various manipulations on the axes ranges so that
        the ranges look OK, e.g., they are 10 percent bigger than
        all graphics object.

        INPUT:
            xmin, xmax, ymin, ymax -- floats that give the axes ranges;
                any can be None, in which case the default to the predefined
                axes ranges for this object.

        OUTPUT:
            good axes ranges

        EXAMPLES:
            sage: P = line([(-1,-2), (3,5)])
            sage: P._prepare_axes(-1,-2, 3, 10)
            (-2.1000000000000001, -0.89000000000000001, 2.2999999999999998, 10.77)
            sage: P._prepare_axes(-1,-2, None, 10)
            (-2.1000000000000001, -0.89000000000000001, -3.2000000000000002, 11.32)
        """
        if xmin is None:
            xmin = self.__xmin
        if xmax is None:
            xmax = self.__xmax
        if ymin is None:
            ymin = self.__ymin
        if ymax is None:
            ymax = self.__ymax

        if xmax < xmin:
            xmax, xmin = xmin, xmax
        elif xmax == xmin:
            x = xmax
            if x == 0:
                x = 1
            xmax = 2*x
            xmin = 0

        if ymax < ymin:
            ymax, ymin = ymin, ymax
        elif ymax == ymin:
            y = ymax
            if y == 0:
                y = 1
            ymax = 2*y
            ymin = 0

        xmin -= 0.1*(xmax-xmin)
        xmax += 0.1*(xmax-xmin)
        ymin -= 0.1*(ymax-ymin)
        ymax += 0.1*(ymax-ymin)

        return xmin,xmax,ymin,ymax

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
        xmin = self.xmin(xmin); xmax = self.xmax(xmax);
        ymin = self.ymin(ymin); ymax = self.ymax(ymax)

        if xmin == xmax:
            xmin -= 1
            xmax += 1
        if ymin == ymax:
            ymin -= 1
            ymax += 1

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
        for g in self.__objects:
            if isinstance(g, GraphicPrimitive_ContourPlot):
                contour = True
            if isinstance(g, GraphicPrimitive_PlotField):
                plotfield = True
            if isinstance(g, GraphicPrimitive_MatrixPlot):
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
            if frame: #add the frame axes
                xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
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
                xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
                # draw the grid
                sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, False)
                # draw the axes
                xmin, xmax, ymin, ymax = sage_axes.add_xy_axes(subplot, xmin, xmax, ymin, ymax)
                subplot.set_xlim(xmin, xmax)
                subplot.set_ylim(ymin, ymax)

            else: #regular plot with no axes
                xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
                subplot.set_xlim(xmin, xmax)
                subplot.set_ylim(ymin, ymax)
                # draw the grid
                sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, False)

        elif (contour or plotfield): #contour or field plot in self.__objects, so adjust axes accordingly
            xmin, xmax = self.__xmin, self.__xmax
            ymin, ymax = self.__ymin, self.__ymax
            subplot.set_xlim([xmin - 0.05*abs(xmax - xmin), xmax + 0.05*abs(xmax - xmin)])
            subplot.set_ylim([ymin - 0.05*abs(ymax - ymin), ymax + 0.05*abs(ymax - ymin)])
            # draw the grid
            sage_gridlines.add_gridlines(subplot, xmin, xmax, ymin, ymax, True)
            # draw the axes
            if axes: #axes=True unless user specifies axes=False
                sage_axes.add_xy_frame_axes(subplot, xmin, xmax, ymin, ymax)

        else: #we have a 'matrix_plot' in self.__objects, so adjust axes accordingly
            xmin, xmax = self.__xmin, self.__xmax
            ymin, ymax = self.__ymin, self.__ymax
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

################## Graphics Primitives ################

class GraphicPrimitive(SageObject):
    """
    Base class for graphics primitives, e.g., things that knows how to draw
    themselves in 2d.

    EXAMPLES:
    We create an object that derives from GraphicPrimitive:
        sage: P = line([(-1,-2), (3,5)])
        sage: P[0]
        Line defined by 2 points
        sage: type(P[0])
        <class 'sage.plot.plot.GraphicPrimitive_Line'>
    """
    def __init__(self, options):
        """
        Create a base class GraphicsPrimitive.  All this does is
        set the options.

        EXAMPLES:
        We indirectly test this function.
            sage: from sage.plot.plot import GraphicPrimitive
            sage: GraphicPrimitive({})
            Graphics primitive
        """

        self.__options = options

    def _allowed_options(self):
        """
        Return the allowed options for a graphics primitive.

        OUTPUT:
            -- a reference to a dictionary.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive
            sage: GraphicPrimitive({})._allowed_options()
            {}
        """
        return {}

    def plot3d(self, **kwds):
        raise NotImplementedError, "3d plotting not implemented for %s" % type(self)

    def _plot3d_options(self, options=None):
        """
        Translate 2d plot options into 3d plot options.
        """
        if options == None:
            options = self.options()
        options_3d = {}
        if 'rgbcolor' in options:
            options_3d['rgbcolor'] = options['rgbcolor']
            del options['rgbcolor']
        if 'alpha' in options:
            options_3d['opacity'] = options['alpha']
            del options['alpha']
        if len(options) != 0:
            raise NotImplementedError, "Unknown plot3d equivalent for %s" % ", ".join(options.keys())
        return options_3d

    def options(self):
        """
        Return the dictionary of options for this graphics primitive.

        By default this function verifies that the options are all
        valid; if any aren't a verbose message is printed with level 0.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive
            sage: GraphicPrimitive({}).options()
            {}
        """
        O = dict(self.__options)
        global do_verify
        if do_verify:
            A = self._allowed_options()
            t = False
            K = A.keys() + ['xmin', 'xmax', 'ymin', 'ymax', 'axes']
            for k in O.keys():
                if not k in K:
                    do_verify = False
                    sage.misc.misc.verbose("WARNING: Ignoring option '%s'=%s"%(k,O[k]), level=0)
                    t = True
            if t:
                s = "\nThe allowed options for %s are:\n"%self
                K.sort()
                for k in K:
                    if A.has_key(k):
                        s += "    %-15s%-60s\n"%(k,A[k])
                sage.misc.misc.verbose(s, level=0)


        if 'hue' in O:
            t = O['hue']
            if not isinstance(t, (tuple,list)):
                t = [t,1,1]
            O['rgbcolor'] = hue(*t)
            del O['hue']
        return O

    def _repr_(self):
        """
        String representation of this graphics primitive.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive
            sage: GraphicPrimitive({})._repr_()
            'Graphics primitive'
        """
        return "Graphics primitive"


class GraphicPrimitive_Arrow(GraphicPrimitive):
    """
    Primitive class that initializes the arrow graphics type

    EXAMPLES:
    We create an arrow graphics object, then take the 0th entry
    in it to get the actual Arrow graphics primitive:
        sage: P = arrow((0,1), (2,3))[0]
        sage: type(P)
        <class 'sage.plot.plot.GraphicPrimitive_Arrow'>
        sage: P
        Arrow from (0.0,1.0) to (2.0,3.0)
    """
    def __init__(self, xtail, ytail, xhead, yhead, options):
        """
        Create an arrow graphics primitive.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive_Arrow
            sage: GraphicPrimitive_Arrow(0,0,2,3,{})
            Arrow from (0.0,0.0) to (2.0,3.0)
        """
        self.xtail = float(xtail)
        self.xhead = float(xhead)
        self.ytail = float(ytail)
        self.yhead = float(yhead)
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        """
        Return the dictionary of allowed options for the arrow graphics primitive.

        EXAMPLES:
             sage: from sage.plot.plot import GraphicPrimitive_Arrow
             sage: list(sorted(GraphicPrimitive_Arrow(0,0,2,3,{})._allowed_options().iteritems()))
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
            sage: from sage.plot.plot import GraphicPrimitive_Arrow
            sage: GraphicPrimitive_Arrow(0,0,2,3,{})._repr_()
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

#TODO: make bar_chart more general
class GraphicPrimitive_BarChart(GraphicPrimitive):
    """
    Graphics primitive that represents a bar chart.

    EXAMPLES:
        sage: from sage.plot.plot import GraphicPrimitive_BarChart
        sage: g = GraphicPrimitive_BarChart(range(4), [1,3,2,0], {}); g
        BarChart defined by a 4 datalist
        sage: type(g)
        <class 'sage.plot.plot.GraphicPrimitive_BarChart'>
    """
    def __init__(self, ind, datalist, options):
        """
        Initialize a BarChart primitive.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive_BarChart
            sage: GraphicPrimitive_BarChart(range(3), [10,3,5], {'width':0.7})
            BarChart defined by a 3 datalist
        """
        self.datalist = datalist
        self.ind = ind
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        """
        Return the allowed options with descriptions for this graphics primitive.
        This is used in displaying an error message when the user gives an option
        that doesn't make sense.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive_BarChart
            sage: g = GraphicPrimitive_BarChart(range(4), [1,3,2,0], {})
            sage: list(sorted(g._allowed_options().iteritems()))
            [('hue', 'The color given as a hue.'),
             ('rgbcolor', 'The color as an rgb tuple.'),
             ('width', 'The width of the bars')]
        """
        return {'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'width':'The width of the bars'}

    def _repr_(self):
        """
        Return text representation of this bar chart graphics primitive.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive_BarChart
            sage: g = GraphicPrimitive_BarChart(range(4), [1,3,2,0], {})
            sage: g._repr_()
            'BarChart defined by a 4 datalist'
        """
        return "BarChart defined by a %s datalist"%(len(self.datalist))

    def _render_on_subplot(self, subplot):
        """
        Render this bar chart graphics primitive on a matplotlib subplot object.

        EXAMPLES:
        This rendering happens implicitly when the following command
        is executed:
            sage: bar_chart([1,2,10])
        """
        options = self.options()
        color = options['rgbcolor']
        width = float(options['width'])
        # it is critical to make numpy arrays of type float below,
        # or bar will go boom:
        import numpy
        ind = numpy.array(self.ind, dtype=float)
        datalist = numpy.array(self.datalist, dtype=float)
        subplot.bar(ind, datalist, color=color, width=width)


class GraphicPrimitive_Line(GraphicPrimitive):
    """
    Primitive class that initializes the line graphics type.

    EXAMPLES:
        sage: from sage.plot.plot import GraphicPrimitive_Line
        sage: GraphicPrimitive_Line([1,2,7], [1,5,-1], {})
        Line defined by 3 points
    """
    def __init__(self, xdata, ydata, options):
        """
        Initialize a line graphics primitive.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive_Line
            sage: GraphicPrimitive_Line([-1,2], [17,4], {'thickness':2})
            Line defined by 2 points
        """
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        """
        Displayed the list of allowed line options.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive_Line
            sage: list(sorted(GraphicPrimitive_Line([-1,2], [17,4], {})._allowed_options().iteritems()))
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
             ('thickness', 'How thick the line is.')]
        """
        return {'alpha':'How transparent the line is.',
                'thickness':'How thick the line is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'linestyle':"The style of the line, which is one of '--' (dashed), '-.' (dash dot), '-' (solid), 'steps', ':' (dotted).",
                'marker':"'0' (tickleft), '1' (tickright), '2' (tickup), '3' (tickdown), '' (nothing), ' ' (nothing), '+' (plus), ',' (pixel), '.' (point), '1' (tri_down), '3' (tri_left), '2' (tri_up), '4' (tri_right), '<' (triangle_left), '>' (triangle_right), 'None' (nothing), 'D' (diamond), 'H' (hexagon2), '_' (hline), '^' (triangle_up), 'd' (thin_diamond), 'h' (hexagon1), 'o' (circle), 'p' (pentagon), 's' (square), 'v' (triangle_down), 'x' (x), '|' (vline)",
                'markersize':'the size of the marker in points',
                'markeredgecolor':'the markerfacecolor can be any color arg',
                'markeredgewidth':'the size of the markter edge in points'
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
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        EXAMPLES:
            sage: EllipticCurve('37a').plot(thickness=5).plot3d()
        """
        from sage.plot.plot3d.shapes2 import line3d
        options = self._plot3d_options()
        options.update(kwds)
        return line3d([(x, y, 0) for x, y in zip(self.xdata, self.ydata)], **options)

    def _repr_(self):
        """
        String representation of a line primitive.

        EXAMPLES:
            sage: from sage.plot.plot import GraphicPrimitive_Line
            sage: GraphicPrimitive_Line([-1,2,3,3], [17,4,0,2], {})._repr_()
            'Line defined by 4 points'
        """
        return "Line defined by %s points"%len(self)

    def __getitem__(self, i):
        """
        Extract the i-th element of the line (which is stored as a list of points).

        INPUT:
            i -- an integer between 0 and the number of points minus 1

        OUTPUT:
            a 2-tuple of floats

        EXAMPLES:
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
            i -- an integer between 0 and the number of points on the
                 line minus 1
            point -- a 2-tuple of floats

        EXAMPLES:
        We create a line graphics object $L$ and get ahold of the
        corresponding line graphics primitive.
            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: line_primitive = L[0]; line_primitive
            Line defined by 4 points

        We then set the 0th point to (0,0) instead of (1,2).
            sage: line_primitive[0] = (0,0)
            sage: line_primitive[0]
            (0.0, 0.0)

        Plotting we visibly see the change -- now the line starts at (0,0).
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
        its length:
            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: len(L[0])
            4
        """
        return len(self.xdata)

    def _render_on_subplot(self, subplot):
        """
        Render this line on a matplotlib subplot.

        INPUT:
            subplot -- a matplotlib subplot

        EXAMPLES:
        This implicitly calls this function:
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


class GraphicPrimitive_Circle(GraphicPrimitive):
    """
    Circle graphics primitive.
    """
    def __init__(self, x, y, r, options):
        self.x = x
        self.y = y
        self.r = r
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        return {'alpha':'How transparent the line is.',
                'fill': 'Whether or not to fill the polygon.',
                'thickness':'How thick the border of the polygon is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

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
        n = 50
        dt = float(2*pi/n)
        x, y, r = self.x, self.y, self.r
        xdata = [x+r*cos(t*dt) for t in range(n+1)]
        ydata = [y+r*sin(t*dt) for t in range(n+1)]
        if fill:
            return GraphicPrimitive_Polygon(xdata, ydata, options).plot3d()
        else:
            return GraphicPrimitive_Line(xdata, ydata, options).plot3d()

class GraphicPrimitive_ContourPlot(GraphicPrimitive):
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

    def _allowed_options(self):
        return {'plot_points':'How many points to use for plotting precision',
                'cmap':"""The colormap, one of (autumn, bone, cool, copper,
                       gray, hot, hsv, jet, pink, prism, spring, summer, winter)""",
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
        cmap = options['cmap']
        contours = options['contours']
        #cm is the matplotlib color map module
        from matplotlib import cm
        try:
            cmap = cm.__dict__[cmap]
        except KeyError:
            from matplotlib.colors import LinearSegmentedColormap as C
            possibilities = ', '.join([str(x) for x in cm.__dict__.keys() if \
                                       isinstance(cm.__dict__[x], C)])
            sage.misc.misc.verbose("The possible color maps include: %s"%possibilities, level = 0)
            raise RuntimeError, "Color map %s not known"%cmap

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


class GraphicPrimitive_MatrixPlot(GraphicPrimitive):
    """
    Primitive class that initializes the
    matrix_plot graphics type
    """
    def __init__(self, xy_data_array, xrange, yrange, options):
        self.xrange = xrange
        self.yrange = yrange
        self.xy_data_array = xy_data_array
        self.xy_array_row = len(xy_data_array)
        self.xy_array_col = len(xy_data_array[0])
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        return {'cmap':"""The colormap, one of (autumn, bone, cool, copper,
                gray, hot, hsv, jet, pink, prism, spring, summer, winter)"""}

    def _repr_(self):
        return "MatrixPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        options = self.options()
        cmap = options['cmap']
        #cm is the matplotlib color map module
        from matplotlib import cm
        try:
            cmap = cm.__dict__[cmap]
        except KeyError:
            from matplotlib.colors import LinearSegmentedColormap as C
            possibilities = ', '.join([str(x) for x in cm.__dict__.keys() if \
                                       isinstance(cm.__dict__[x], C)])
            sage.misc.misc.verbose("The possible color maps include: %s"%possibilities, level=0)
            raise RuntimeError, "Color map %s not known"%cmap

        subplot.imshow(self.xy_data_array, cmap=cmap, interpolation='nearest', extent=(0,self.xrange[1],0,self.yrange[1]))

# Below is the base class that is used to make 'field plots'.
# Its implementation is motivated by 'PlotField'.
# Currently it is used to make the function 'plot_vector_field'
# TODO: use this to make these functions:
# 'plot_gradient_field' and 'plot_hamiltonian_field'

class GraphicPrimitive_PlotField(GraphicPrimitive):
    """
    Primitive class that initializes the
    plot_field graphics type
    """
    def __init__(self, xpos_array, ypos_array, xvec_array, yvec_array, options):
        self.xpos_array = xpos_array
        self.ypos_array = ypos_array
        self.xvec_array = xvec_array
        self.yvec_array = yvec_array
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        return {'plot_points':'How many points to use for plotting precision',
                'pivot': 'Where the arrow should be placed in relation to the point (tail, middle, tip)',
                'headwidth': 'Head width as multiple of shaft width, default is 3',
                'headlength': 'head length as multiple of shaft width, default is 5',
                'headaxislength': 'head length at shaft intersection, default is 4.5'}

    def _repr_(self):
        return "PlotField defined by a %s x %s vector grid"%(len(self.xpos_array), len(self.ypos_array))

    def _render_on_subplot(self, subplot):
        options = self.options()
        quiver_options = options.copy()
        quiver_options.pop('plot_points')
        subplot.quiver(self.xpos_array, self.ypos_array, self.xvec_array, self.yvec_array, **quiver_options)

class GraphicPrimitive_Disk(GraphicPrimitive):
    """
    Primitive class that initializes the
    disk graphics type
    """
    def __init__(self, point, r, angle, options):
        self.x = point[0]
        self.y = point[1]
        self.r = r
        self.rad1 = angle[0]
        self.rad2 = angle[1]
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        return {'alpha':'How transparent the line is.',
                'fill': 'Whether or not to fill the polygon.',
                'thickness':'How thick the border of the polygon is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

    def _repr_(self):
        return "Disk defined by (%s,%s) with r=%s spanning (%s, %s) radians"%(self.x,
        self.y, self.r, self.rad1, self.rad2)

    def _render_on_subplot(self, subplot):
        import matplotlib.patches as patches
        options = self.options()
        deg1 = self.rad1*(360.0/(2.0*pi)) #convert radians to degrees
        deg2 = self.rad2*(360.0/(2.0*pi))
        p = patches.Wedge((float(self.x), float(self.y)), float(self.r), float(deg1),
                            float(deg2))
        p.set_linewidth(float(options['thickness']))
        p.set_fill(options['fill'])
        p.set_alpha(options['alpha'])
        c = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(c)
        p.set_facecolor(c)
        subplot.add_patch(p)

class GraphicPrimitive_Point(GraphicPrimitive):
    """
    Primitive class that initializes the
    point graphics type

    """
    def __init__(self, xdata, ydata, options):
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        return {'alpha':'How transparent the line is.',
                'pointsize': 'How big the point is.',
                'faceted': 'If True color the edge of the point.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

    def _plot3d_options(self, options=None):
        if options == None:
            options = dict(self.options())
        options_3d = {}
        if 'pointsize' in options:
            options_3d['size'] = options['pointsize']
            del options['pointsize']
        if 'faceted' in options:
            if options['faceted']:
                raise NotImplementedError, "No 3d faceted points."
            del options['faceted']
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: P = E(0,0)
            sage: def get_points(n):
            ...     return sum([point((i*P)[0:2], pointsize=3) for i in range(-n,n) if i != 0 and (i*P)[0] < 3])
            sage: sum([get_points(15*n).plot3d(z=n) for n in range(1,10)])
        """
        from sage.plot.plot3d.base import Graphics3dGroup
        from sage.plot.plot3d.shapes2 import point3d
        options = self._plot3d_options()
        options.update(kwds)
        all = [point3d([(x, y, 0) for x, y in zip(self.xdata, self.ydata)], **options)]
        if len(all) == 1:
            return all[0]
        else:
            return Graphics3dGroup(all)

    def _repr_(self):
        return "Point set defined by %s point(s)"%len(self.xdata)

    def __getitem__(self, i):
        return self.xdata[i], self.ydata[i]

    def _render_on_subplot(self,subplot):
        """
        TESTS:
        We check to make sure that #2076 is fixed by verifying all
        the points are red.
            sage: point(((1,1), (2,2), (3,3)), rgbcolor=hue(1), pointsize=30)
        """
        options = self.options()

        #Convert the color to a hex string so that the scatter
        #method does not interpret it as a list of 3 floating
        #point color specifications when there are
        #three points. This is mentioned in the matplotlib 0.98
        #documentation and fixes #2076
        from matplotlib.colors import rgb2hex
        c = rgb2hex(to_mpl_color(options['rgbcolor']))

        a = float(options['alpha'])
        s = int(options['pointsize'])
        faceted = options['faceted'] #faceted=True colors the edge of point
        scatteroptions={}
        if not faceted: scatteroptions['edgecolors'] = 'none'

        subplot.scatter(self.xdata, self.ydata, s=s, c=c, alpha=a, **scatteroptions)

class GraphicPrimitive_Polygon(GraphicPrimitive):
    """
    Primitive class that initializes the
    polygon graphics type

    """
    def __init__(self, xdata, ydata, options):
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive.__init__(self, options)

    def _repr_(self):
        return "Polygon defined by %s points"%len(self)

    def __getitem__(self, i):
        return self.xdata[i], self.ydata[i]

    def __setitem__(self, i, point):
        i = int(i)
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        return len(self.xdata)

    def _allowed_options(self):
        return {'alpha':'How transparent the line is.',
                'thickness': 'How thick the border line is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

    def _plot3d_options(self, options=None):
        if options == None:
            options = dict(self.options())
        if 'thickness' in options:
            del options['thickness']
        return GraphicPrimitive._plot3d_options(self, options)

    def plot3d(self, **kwds):
        """
        EXAMPLES:
            sage: polygon([(cos(t), sin(t)) for t in srange(0, 2*pi, 2*pi/5)]).plot3d()
        """
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        options = self._plot3d_options()
        options.update(kwds)
        return IndexFaceSet([[(x, y, 0) for x, y in zip(self.xdata, self.ydata)]], **options)

    def _render_on_subplot(self, subplot):
        import matplotlib.patches as patches
        options = self.options()
        p = patches.Polygon([(self.xdata[i],self.ydata[i]) for i in xrange(len(self.xdata))])
        p.set_linewidth(float(options['thickness']))
        a = float(options['alpha'])
        p.set_alpha(a)
        c = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(c)
        p.set_facecolor(c)
        subplot.add_patch(p)

class GraphicPrimitive_Text(GraphicPrimitive):
    """
    Text graphics primitive.
    """
    def __init__(self, string, point, options):
        self.string = string
        self.x = point[0]
        self.y = point[1]
        GraphicPrimitive.__init__(self, options)

    def _repr_(self):
        return "Text %s at the point (%s,%s)"%(self.string, self.x, self.y)

    def _allowed_options(self):
        return {'fontsize': 'How big the text is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'axis_coords':'Uses axis coordinates -- (0,0) lower left and (1,1) upper right',
                'vertical_alignment': 'how to align vertically: top, center, bottom',
                'horizontal_alignment':'how to align horizontally: left, center, right'}

    def _plot3d_options(self, options=None):
        if options == None:
            options = dict(self.options())
        options_3d = {}
        # TODO: figure out how to implement rather than ignore
        if 'fontsize' in options:
            del options['fontsize']
        if 'vertical_alignment' in options:
            del options['vertical_alignment']
        if 'horizontal_alignment' in options:
            del options['horizontal_alignment']
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        from sage.plot.plot3d.shapes2 import text3d
        options = self._plot3d_options()
        options.update(kwds)
        return text3d(self.string, (self.x, self.y, 0), **options)

    def _render_on_subplot(self, subplot):
        options = self.options()
        opts = {}
        opts['color'] = options['rgbcolor']
        opts['fontsize'] = int(options['fontsize'])
        opts['verticalalignment'] = options['vertical_alignment']
        opts['horizontalalignment'] = options['horizontal_alignment']
        if options['axis_coords']:
            opts['transform'] = subplot.transAxes
        subplot.text(float(self.x), float(self.y), self.string, **opts)

class GraphicPrimitive_NetworkXGraph(GraphicPrimitive):
    """
    Primitive class that initializes the NetworkX graph type.

    INPUT:
        graph -- a NetworkX graph
        pos -- an optional positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {   0: [-0.91679746, 0.88169588,],
                1: [ 0.47294849, 1.125     ,],
                2: [ 1.125     ,-0.12867615,],
                3: [ 0.12743933,-1.125     ,],
                4: [-1.125     ,-0.50118505,]   }
        vertex_labels -- determines whether labels for nodes are plotted
        vertex_size -- node size
        vertex_colors -- a dictionary specifying node colors: each key is a color recognized by
                        matplotlib, and each entry is a list of vertices.
        edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                        matplotlib, and each entry is a list of edges.
        scaling_term -- default is 0.05. if nodes are getting chopped off, increase; if graph
                        is too small, decrease. should be positive, but values much bigger than
                        1/8 won't be useful unless the nodes are huge
        draw_edges -- whether to draw edges.

    EXAMPLES:
        sage: from sage.plot.plot import GraphicPrimitive_NetworkXGraph
        sage: import networkx
        sage: D = networkx.dodecahedral_graph()
        sage: NGP = GraphicPrimitive_NetworkXGraph(D)
        sage: g = Graphics()
        sage: g._Graphics__objects.append(NGP)
        sage: g.axes(False)
        sage: g.show()

        sage: import networkx
        sage: from sage.plot.plot import GraphicPrimitive_NetworkXGraph
        sage: from math import sin, cos, pi
        sage: P = networkx.petersen_graph()
        sage: d = {'#FF0000':[0,5], '#FF9900':[1,6], '#FFFF00':[2,7], '#00FF00':[3,8], '#0000FF':[4,9]}
        sage: pos_dict = {}
        sage: for i in range(5):
        ...    x = float(cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: for i in range(10)[5:]:
        ...    x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: NGP = GraphicPrimitive_NetworkXGraph(graph=P, vertex_colors=d, pos=pos_dict)
        sage: g = Graphics()
        sage: g._Graphics__objects.append(NGP)
        sage: g.axes(False)
        sage: g.show()

        sage: from sage.plot.plot import rainbow
        sage: from sage.plot.plot import GraphicPrimitive_NetworkXGraph
        sage: import networkx
        sage: C = graphs.CubeGraph(5)
        sage: pos = C.get_pos()
        sage: G = C.networkx_graph()
        sage: R = rainbow(5)
        sage: edge_colors = {}
        sage: for i in range(5):
        ...    edge_colors[R[i]] = []
        sage: for u,v,l in C.edges():
        ...    for i in range(5):
        ...        if u[i] != v[i]:
        ...            edge_colors[R[i]].append((u,v,l))
        sage: NGP = GraphicPrimitive_NetworkXGraph(G, pos=pos, vertex_labels=False, vertex_size=0, edge_colors=edge_colors)
        sage: G = Graphics()
        sage: G._Graphics__objects.append(NGP)
        sage: G.axes_range(xmin=-1.1, xmax=2.2, ymin=0, ymax=3.25)
        sage: G.axes(False)
        sage: G.show()

    We color the edges and vertices of a Dodecahedral graph:
        sage: g = graphs.DodecahedralGraph()
        sage: g.show(edge_colors={(1.0, 0.8571428571428571, 0.0): g.edges()})

    """
    def __init__(self, graph, pos=None, vertex_labels=True, vertex_size=300, \
                   vertex_colors=None, edge_colors=None, scaling_term=0.05, \
                   draw_edges=True):
        self.__nxg = graph
        self.__vertex_size = vertex_size
        self.__vertex_labels = vertex_labels
        self.__vertex_colors = vertex_colors
        self.__edge_colors = edge_colors
        self.__draw_edges = draw_edges
        if len(self.__nxg) != 0:
            import networkx as NX
            if pos is None:
                self.__pos = NX.drawing.spring_layout(self.__nxg)
            else:
                self.__pos = pos

            nodelist=self.__nxg.nodes()

            xes = [self.__pos[v][0] for v in self.__pos]
            ys = [self.__pos[v][1] for v in self.__pos]
            xmin = min(xes)
            xmax = max(xes)
            ymin = min(ys)
            ymax = max(ys)
            if xmax == xmin:
                xmax += 1
                xmin -= 1
            if ymax == ymin:
                ymax += 1
                ymin -= 1
            dx = xmax - xmin
            dy = ymax - ymin

            if not pos is None:
                random = current_randstate().python_random().random
                missing = []
                for v in nodelist:
                    if not v in self.__pos:
                        missing.append(v)
                for v in missing:
                    self.__pos[v] = [xmin + dx*random(),ymin + dy*random()]

            # adjust the plot
            xmin -= scaling_term*dx
            xmax += scaling_term*dx
            ymin -= scaling_term*dy
            ymax += scaling_term*dy

            self._xmin = xmin
            self._xmax = xmax
            self._ymin = ymin
            self._ymax = ymax
        else:
            self.__pos = {}
            self._xmin = -1
            self._xmax = 1
            self._ymin = -1
            self._ymax = 1

    def _render_on_subplot(self, subplot):
        if len(self.__nxg) != 0:
            import networkx as NX
            vertex_size = float(self.__vertex_size)
            if self.__vertex_colors is None:
                NX.draw_networkx_nodes(G=self.__nxg, pos=self.__pos, ax=subplot, node_size=vertex_size)
            else:
                from numpy import array
                for i in self.__vertex_colors:
                    NX.draw_networkx_nodes(G=self.__nxg, nodelist=self.__vertex_colors[i],
                                           node_color=i if isinstance(i, str) else [float(z) for z in i],
                                           pos=self.__pos, ax=subplot, node_size=vertex_size)
            if self.__draw_edges:
                if self.__edge_colors is None:
                    NX.draw_networkx_edges(G=self.__nxg, pos=self.__pos, ax=subplot, node_size=vertex_size)
                else:
                    for i in self.__edge_colors:
                        NX.draw_networkx_edges(G=self.__nxg, pos=self.__pos, edgelist=self.__edge_colors[i],
                                               edge_color=i if isinstance(i, str) else [float(z) for z in i],
                                               ax=subplot, node_size=vertex_size)
            if self.__vertex_labels:
                labels = {}
                for v in self.__nxg:
                    labels[v] = str(v)
                NX.draw_networkx_labels(self.__nxg, self.__pos, labels=labels, ax=subplot)



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

def arrow(tailpoint, headpoint, **kwds):
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


    TESTS:
    We check to make sure the x/y min/max data is correct:
        sage: a = arrow((1,1), (5,5))
        sage: a.xmin()
        1.0
        sage: a.xmax()
        5.0
    """

    options = {'width':2,'rgbcolor':(0, 0, 1)}
    options.update(kwds)

    #Get the x/y min/max data
    xtail = float(tailpoint[0])
    ytail = float(tailpoint[1])
    xhead = float(headpoint[0])
    yhead = float(headpoint[1])
    xmin = min(xtail, xhead)
    xmax = max(xtail, xhead)
    ymin = min(ytail, yhead)
    ymax = max(ytail, yhead)

    g = Graphics(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    g._arrow(xtail, ytail, xhead, yhead, options=options)
    return g

def bar_chart(datalist, **kwds):
    """
    A bar chart of (currently) one list of numerical data.
    Support for more datalists in progress.

    EXAMPLES:
    A bar_chart with blue bars:
        sage: bar_chart([1,2,3,4])

    A bar_chart with thinner bars:
        sage: bar_chart([x^2 for x in range(1,20)], width=0.2)

    A bar_chart with negative values and red bars:
        sage: bar_chart([-3,5,-6,11], rgbcolor=(1,0,0))
    """
    options = {'width':0.5,'rgbcolor':(0, 0, 1)}
    options.update(kwds)

    dl = len(datalist)
    #if dl > 1:
    #    print "WARNING, currently only 1 data set allowed"
    #    datalist = datalist[0]
    if dl == 3:
        datalist = datalist+[0]
    #bardata = []
    #cnt = 1
    #for pnts in datalist:
        #ind = [i+cnt/dl for i in range(len(pnts))]
    ind = range(len(datalist))
    xrange = (0, len(datalist))
    yrange = (min(datalist), max(datalist))
        #bardata.append([ind, pnts, xrange, yrange])
        #cnt += 1

    g = Graphics()
    #TODO: improve below for multiple data sets!
    #cnt = 1
    #for ind, pnts, xrange, yrange in bardata:
        #options={'rgbcolor':hue(cnt/dl),'width':0.5/dl}
    #    g._bar_chart(ind, pnts, xrange, yrange, options=options)
    #    cnt += 1
    #else:
    g._bar_chart(ind, datalist, xrange, yrange, options=options)
    return g


def circle(point, radius, **kwds):
    """
    Return a circle at a point = $(x,y)$ with radius = $r$.

    circle(center, radius, **kwds)

    INPUT:
        center -- a 2-tuple (x,y)
        radius -- a positive number
        alpha -- default: 1
        fill -- default: False
        thickness -- default: 1
        rgbcolor -- default: (0,0,0)

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

    TESTS:
    We test to make sure the x/y min/max data is set correctly.
        sage: p = circle((3, 3), 1)
        sage: p.xmin()
        2.0
        sage: p.ymin()
        2.0
    """
    options={'alpha':1,'fill':False,'thickness':1,'rgbcolor':(0, 0, 1)}
    options.update(kwds)

    r = float(radius)
    point = (float(point[0]), float(point[1]))
    g = Graphics(xmin=point[0]-r, xmax=point[0]+r, ymin=point[1]-r, ymax=point[1]+r)
    g._circle(point[0], point[1], r, options)
    return g

def contour_plot(f, xrange, yrange, **kwds):
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
        cmap         -- string (default: 'gray'), the color map to use:
                        autumn, bone, cool, copper, gray, hot, hsv,
                        jet, pink, prism, spring, summer, winter
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
        sage: contour_plot(f, (-2, 2), (-2, 2), contours=2)
        sage: contour_plot(f, (-2, 2), (-2, 2), contours=(0.1, 1.0, 1.2, 1.4), cmap='hsv')
        sage: contour_plot(f, (-2, 2), (-2, 2), contours=(1.0,), fill=False)


    TESTS:
    We test to make sure that the x/y min/max data is set correctly.
        sage: p = contour_plot(f, (3, 6), (3, 6))
        sage: p.xmin()
        3.0
        sage: p.ymin()
        3.0
    """
    options = {'plot_points':25, 'fill':True, 'cmap':'gray', 'contours':None}
    options.update(kwds)

    g, xstep, ystep, xrange, yrange = setup_for_eval_on_grid([f], xrange, yrange, options['plot_points'])
    g = g[0]
    xy_data_array = [[g(x, y) for x in \
                      sage.misc.misc.xsrange(xrange[0], xrange[1], xstep)]
                      for y in sage.misc.misc.xsrange(yrange[0], yrange[1], ystep)]

    g = Graphics(xmin=float(xrange[0]), xmax=float(xrange[1]), ymin=float(yrange[0]), ymax=float(yrange[1]))
    g._contour_plot(xy_data_array, xrange, yrange, options)
    return g

def implicit_plot(f, xrange, yrange, **kwds):
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


    EXAMPLES:

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
    options = {'plot_points':25, 'fill':False, 'cmap':'gray', 'contours':(0.0,)}
    options.update(kwds)
    return contour_plot(f, xrange, yrange, **options)

def line(points, **kwds):
    """
    Returns either a 2-dimensional or 3-dimensional line depending
    on value of points.

    For information regarding additional arguments, see either line2d?
    or line3d?.

    EXAMPLES:
        sage: line([(0,0), (1,1)])
        sage: line([(0,0,1), (1,1,1)])
    """
    try:
        return line2d(points, **kwds)
    except ValueError:
        from sage.plot.plot3d.shapes2 import line3d
        return line3d(points, **kwds)

def line2d(points, **kwds):
    r"""
    Create the line through the given list of points.

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
       markeredgewidth -- the size of the markter edge in points

    EXAMPLES:
    A blue conchoid of Nicomedes:

        sage: L = [[1+5*cos(pi/2+pi*i/100), tan(pi/2+pi*i/100)*(1+5*cos(pi/2+pi*i/100))] for i in range(1,100)]
        sage: line(L, rgbcolor=(1/4,1/8,3/4))

    A line with 2 complex points:
        sage: i = CC.0
        sage: line([1+i, 2+3*i])

    A blue hypotrochoid (3 leaves):
        sage: n = 4; h = 3; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: line(L, rgbcolor=(1/4,1/4,3/4))

    A blue hypotrochoid (4 leaves):

        sage: n = 6; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: line(L, rgbcolor=(1/4,1/4,3/4))

    A red limacon of Pascal:

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,101)]
        sage: line(L, rgbcolor=(1,1/4,1/2))

    A light green trisectrix of Maclaurin:

        sage: L = [[2*(1-4*cos(-pi/2+pi*i/100)^2),10*tan(-pi/2+pi*i/100)*(1-4*cos(-pi/2+pi*i/100)^2)] for i in range(1,100)]
        sage: line(L, rgbcolor=(1/4,1,1/8))

    A green lemniscate of Bernoulli:

        sage: v = [(1/cos(-pi/2+pi*i/100), tan(-pi/2+pi*i/100)) for i in range(201)]
        sage: L = [(a/(a^2+b^2), b/(a^2+b^2)) for a,b in v]
        sage: line(L, rgbcolor=(1/4,3/4,1/8))

    A red plot of the Jacobi elliptic function $\text{sn}(x,2)$, $-3<x<3$:

        sage: L = [(i/100.0, jacobi('sn', i/100.0 ,2.0)) for i in range(-300,300,30)]
        sage: line(L, rgbcolor=(3/4,1/4,1/8))

    A red plot of $J$-Bessel function $J_2(x)$, $0<x<10$:

        sage: L = [(i/10.0, bessel_J(2,i/10.0)) for i in range(100)]
        sage: line(L, rgbcolor=(3/4,1/4,5/8))


    A purple plot of the Riemann zeta function $\zeta(1/2 + it)$, $0<t<30$:

        sage: i = CDF.gen()
        sage: v = [zeta(0.5 + n/10 * i) for n in range(300)]
        sage: L = [(z.real(), z.imag()) for z in v]
        sage: line(L, rgbcolor=(3/4,1/2,5/8))

    A purple plot of the Hasse-Weil $L$-function $L(E, 1 + it)$, $-1<t<10$:

        sage: E = EllipticCurve('37a')
        sage: vals = E.lseries().values_along_line(1-I, 1+10*I, 100) # critical line
        sage: L = [(z[1].real(), z[1].imag()) for z in vals]
        sage: line(L, rgbcolor=(3/4,1/2,5/8))

    A red, blue, and green "cool cat":

        sage: G = plot(-cos(x), -2, 2, thickness=5, rgbcolor=(0.5,1,0.5))
        sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))
        sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
        sage: G + P + Q   # show the plot

    A line with no points or one point:
        sage: line([])
        sage: line([(1,1)])

    TESTS:
    We test to make sure that the x/y min/max data are set correctly.
        sage: l = line([(100, 100), (120, 120)])
        sage: l.xmin()
        100.0
        sage: l.xmax()
        120.0

    """
    options = {'alpha':1,'rgbcolor':(0,0,1),'thickness':1}
    options.update(kwds)

    xdata, ydata = xydata_from_point_list(points)
    g = Graphics(**minmax_data(xdata, ydata, dict=True))
    g._Graphics__objects.append(GraphicPrimitive_Line(xdata, ydata, options))
    return g


def matrix_plot(mat, **kwds):
    r"""
    A plot of a given matrix or 2D array.

    Each ($i$th, $j$th) matrix element is given a different
    color value depending on its relative size compared
    to the other elements in the matrix.

    The tick marks drawn on the frame axes denote the
    ($i$th, $j$th) element of the matrix.

    EXAMPLES:

    A matrix over ZZ colored with different grey levels:

        sage: matrix_plot(matrix([[1,3,5,1],[2,4,5,6],[1,3,5,7]]))

    Here we make a random matrix over RR and use cmap='hsv'
    to color the matrix elements different RGB colors:

        sage: matrix_plot(random_matrix(RDF, 50), cmap='hsv')

    Another random plot, but over GF(389):
        sage: matrix_plot(random_matrix(GF(389), 10), cmap='Oranges')

    """
    options = {'cmap':'gray'}
    options.update(kwds)

    from sage.matrix.all import is_Matrix
    from matplotlib.numerix import array
    if not is_Matrix(mat) or (isinstance(mat, (list, tuple)) and isinstance(mat[0], (list, tuple))):
        raise TypeError, "mat must be of type Matrix or a two dimensional array"

    if is_Matrix(mat):
        xrange = (0, mat.ncols())
        yrange = (0, mat.nrows())
    else:
        xrange = (0, len(mat[0]))
        yrange = (0, len(mat))
    xy_data_array = [array(r, dtype=float) for r in mat]

    g = Graphics()
    g._matrix_plot(xy_data_array, xrange, yrange, options)
    return g




# Below is the base class that is used to make 'plot_vector_field'.
# Its implementation is motivated by 'PlotVectorField'.
# TODO: make class similiar to this one to implement:
# 'plot_gradient_field' and 'plot_hamiltonian_field'
def plot_vector_field((f, g), xrange, yrange, **kwds):
    r"""

    \code{plot_vector_field} takes two functions of two variables, $(f(x,y), g(x,y))$
    and plots vector arrows of the function over the specified
    xrange and yrange as demonstrated below.

    plot_vector_field((f, g), (xvar, xmin, xmax), (yvar, ymin, ymax))

    EXAMPLES:
    Plot the vector fields involving sin and cos
        sage: x,y = var('x y')
        sage: plot_vector_field((sin(x), cos(y)), (x,-3,3), (y,-3,3))
        sage: plot_vector_field(( y, (cos(x)-2)*sin(x)), (x,-pi,pi), (y,-pi,pi))

    Plot a gradient field
        sage: u,v = var('u v')
        sage: f = exp(-(u^2+v^2))
        sage: plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2))


    TESTS:
        sage: p = plot_vector_field((.01*x,x+y), (10,20), (10,20))
        sage: p.xmin()
        10.0
        sage: p.ymin()
        10.0

    """
    options = {'plot_points':20}
    options.update(kwds)

    z, xstep, ystep, xrange, yrange = setup_for_eval_on_grid([f,g], xrange, yrange, options['plot_points'])
    f,g = z

    xpos_array, ypos_array, xvec_array, yvec_array = [],[],[],[]
    for x in sage.misc.misc.xsrange(xrange[0], xrange[1], xstep):
        for y in sage.misc.misc.xsrange(yrange[0], yrange[1], ystep):
            xpos_array.append(x)
            ypos_array.append(y)
            xvec_array.append(f(x,y))
            yvec_array.append(g(x,y))

    import numpy
    xvec_array = numpy.array(xvec_array, dtype=float)
    yvec_array = numpy.array(yvec_array, dtype=float)
    g = Graphics(xmin=xrange[0], xmax=xrange[1], ymin=yrange[0],  ymax=yrange[1])
    g._plot_field(xpos_array, ypos_array, xvec_array, yvec_array, xrange, yrange, options)
    return g

def plot_slope_field(f, xrange, yrange, **kwds):
    r"""

    \code{plot_slope_field} takes a function of two variables, $f(x,y)$, and at various points (x_i,y_i), plots a line with slope $f(x_i,y_i)$

    plot_slope_field((f, g), (xvar, xmin, xmax), (yvar, ymin, ymax))

    EXAMPLES:
    A logistic function modeling population growth.
        sage: x,y = var('x y')
        sage: capacity = 3 # thousand
        sage: growth_rate = 0.7 # population increases by 70% per unit of time
        sage: plot_slope_field(growth_rate*(1-y/capacity)*y, (x,0,5), (y,0,capacity*2)).show(aspect_ratio=1)

    Plot a slope field involving sin and cos
        sage: x,y = var('x y')
        sage: plot_slope_field(sin(x+y)+cos(x+y), (x,-3,3), (y,-3,3)).show(aspect_ratio=1)

    """
    from math import sqrt
    slope_options = {'headaxislength': 0, 'headlength': 0, 'pivot': 'middle'}
    slope_options.update(kwds)

    from sage.calculus.calculus import sqrt
    norm = sqrt((f**2+1))
    return plot_vector_field((1/norm, f/norm), xrange, yrange, **slope_options)

def disk(point, radius, angle, **kwds):
    r"""

    A disk at a point = $(x,y)$ with radius = $r$
    spanning (in radians) angle=$(rad1, rad2)$

    EXAMPLES:
    Make some dangerous disks:

        sage: bl = disk((0.0,0.0), 1, (pi, 3*pi/2), rgbcolor=(1,1,0))
        sage: tr = disk((0.0,0.0), 1, (0, pi/2), rgbcolor=(1,1,0))
        sage: tl = disk((0.0,0.0), 1, (pi/2, pi), rgbcolor=(0,0,0))
        sage: br = disk((0.0,0.0), 1, (3*pi/2, 2*pi), rgbcolor=(0,0,0))
        sage: P  = tl+tr+bl+br
        sage: P.show(figsize=(4,4),xmin=-2,xmax=2,ymin=-2,ymax=2)

    TESTS:
    Check to make sure that the x/y min/max data is correctly set.
        sage: d = disk((5,5), 1, (pi/2, pi), rgbcolor=(0,0,0))
        sage: d.xmin()
        4.0
        sage: d.ymin()
        4.0
        sage: d.xmax()
        6.0
    """
    options = {'alpha':1,'fill':True,'rgbcolor':(0,0,1),'thickness':0}
    options.update(kwds)

    r = float(radius)
    point = (float(point[0]), float(point[1]))
    angle = (float(angle[0]), float(angle[1]))

    xmin = point[0] - r
    xmax = point[0] + r
    ymin = point[1] - r
    ymax = point[1] + r
    g = Graphics(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    g._disk(point, r, angle, options)
    return g

def point(points, **kwds):
    """
    Returns either a 2-dimensional or 3-dimensional point depending
    on value of points.

    For information regarding additional arguments, see either point2d?
    or point3d?.

    EXAMPLES:
        sage: point([(0,0), (1,1)])
        sage: point([(0,0,1), (1,1,1)])
    """
    try:
        return point2d(points, **kwds)
    except ValueError:
        from sage.plot.plot3d.shapes2 import point3d
        return point3d(points, **kwds)


def point2d(points, **kwds):
    r"""

    A point of size `pointsize' defined by point = $(x,y)$.
    Point takes either a single tuple of coordinates or a list of tuples.

    EXAMPLES:
        A purple point from a single tuple or coordinates:
        sage: point((0.5, 0.5), rgbcolor=hue(0.75))

        Here are some random larger red points, given as a list of tuples
        sage: point(((0.5, 0.5), (1, 2), (0.5, 0.9), (-1, -1)), rgbcolor=hue(1), pointsize=30)

    TESTS:
    We check to make sure that the x/y min/max data is set correctly.
        sage: p = point((3, 3), rgbcolor=hue(0.75))
        sage: p.xmin()
        3.0
        sage: p.ymin()
        3.0

    """
    options = {'alpha':1,'pointsize':10,'faceted':False,'rgbcolor':(0,0,1)}
    options.update(kwds)

    xdata, ydata = xydata_from_point_list(points)
    g = Graphics(**minmax_data(xdata, ydata, dict=True))
    g._point(xdata, ydata, options, extend_axes=False)
    return g

points = point

def polygon(points, **kwds):
    r"""

    EXAMPLES:
    We create a purple-ish polygon:
        sage: polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))

    Some modern art -- a random polygon:
        sage: v = [(randrange(-5,5), randrange(-5,5)) for _ in range(10)]
        sage: polygon(v)

    A purple hexagon:

        sage: L = [[cos(pi*i/3),sin(pi*i/3)] for i in range(6)]
        sage: polygon(L, rgbcolor=(1,0,1))

    A green deltoid:

        sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
        sage: polygon(L, rgbcolor=(1/8,3/4,1/2))

    A blue hypotrochoid:

        sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
        sage: polygon(L, rgbcolor=(1/8,1/4,1/2))

    Another one:

        sage: n = 4; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: polygon(L, rgbcolor=(1/8,1/4,3/4))

    A purple epicycloid:

        sage: m = 9; b = 1
        sage: L = [[m*cos(pi*i/100)+b*cos((m/b)*pi*i/100),m*sin(pi*i/100)-b*sin((m/b)*pi*i/100)] for i in range(200)]
        sage: polygon(L, rgbcolor=(7/8,1/4,3/4))

    A brown astroid:

        sage: L = [[cos(pi*i/100)^3,sin(pi*i/100)^3] for i in range(200)]
        sage: polygon(L, rgbcolor=(3/4,1/4,1/4))

    And, my favorite, a greenish blob:

        sage: L = [[cos(pi*i/100)*(1+cos(pi*i/50)), sin(pi*i/100)*(1+sin(pi*i/50))] for i in range(200)]
        sage: polygon(L, rgbcolor=(1/8, 3/4, 1/2))

    This one is for my wife:

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,100)]
        sage: polygon(L, rgbcolor=(1,1/4,1/2))

    TESTS:
    We check to make sure that the x/y min/max data is set correctly.
        sage: p = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))
        sage: p.ymin()
        0.0
        sage: p.xmin()
        1.0


    AUTHORS:
        -- David Joyner (2006-04-14): the long list of examples above.

    """
    options = {'alpha':1,'rgbcolor':(0,0,1),'thickness':0}
    options.update(kwds)

    xdata, ydata = xydata_from_point_list(points)
    g = Graphics(**minmax_data(xdata, ydata, dict=True))
    g._polygon(xdata, ydata, options, extend_axes=False)
    return g

def plot(funcs, *args, **kwds):
    r"""
    Use plot by writing

        \code{plot(X, ...)}

    where $X$ is a \sage object (or list of \sage objects) that either is
    callable and returns numbers that can be coerced to floats, or has
    a plot method that returns a \class{GraphicPrimitive} object.

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

        sage: p = plot([sin(x), cos(x)], 100, 120)
        sage: p.xmin()
        100.0
        sage: p.xmax()
        120.0

    We check to handle cases where the function gets evaluated at a point
    which causes an 'inf' or '-inf' result to be produced.
        sage: p = plot(1/x, 0, 1)
        sage: p = plot(-1/x, 0, 1)
    """
    do_show = kwds.pop('show',False)
    if hasattr(funcs, 'plot'):
        G = funcs.plot(*args, **kwds)
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
    if do_show:
        G.show()
    return G

def _plot(funcs, xrange, parametric=False,
              polar=False, label='', randomize=True, **kwds):
    options = {'alpha':1,'thickness':1,'rgbcolor':(0,0,1),
               'plot_points':200, 'adaptive_tolerance':0.01,
               'adaptive_recursion':5, 'rgbcolor': (0,0,1) }

    if kwds.has_key('color') and not kwds.has_key('rgbcolor'):
        kwds['rgbcolor'] = kwds['color']
        del kwds['color']

    options.update(kwds)

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
        return reduce(operator.add, (plot(f, (xmin, xmax), polar=polar, **kwds) for f in funcs))

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
    G = line(data, **options)

    # Label?
    if label:
        label = '  '+str(label)
        G += text(label, data[-1], horizontal_alignment='left',
                  vertical_alignment='center')

    return G


def text(string, (x,y), **kwds):
    r"""
    Returns a 2d text graphics object at the point $(x,y)$.

    2D OPTIONS:
        fontsize -- How big the text is
        rgbcolor -- The color as an rgb tuple
        hue -- The color given as a hue
        vertical_alignment -- how to align vertically: top, center, bottom
        horizontal_alignment -- how to align horizontally: left, center, right
        axis_coords -- (default: False) if True, use axis coordinates, so that
                       (0,0) is the lower left and (1,1) upper right, irregardless
                       of the x and y range of plotted values.

    EXAMPLES:
    Some text:
        sage: text("Sage is really neat!!",(2,12))

    The same text in larger font and colored red:
        sage: text("Sage is really neat!!",(2,12),fontsize=20,rgbcolor=(1,0,0))

    Some text but guaranteed to be in the lower left no matter what:
        sage: text("Sage is really neat!!",(0,0), axis_coords=True, horizontal_alignment='left')

    You can also align text differently:
        sage: t1 = text("Hello",(1,1), vertical_alignment="top")
        sage: t2 = text("World", (1,0.5), horizontal_alignment="left")
        sage: t1 + t2   # render the sume
    """
    options = {'fontsize':10, 'rgbcolor':(0,0,1),
               'horizontal_alignment':'center',
               'vertical_alignment':'center',
               'axis_coords':False}
    options.update(kwds)

    point = (float(x), float(y))
    g = Graphics()
    g._text(string, point, options)
    return g


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
        sage: p = list_plot([(100,100), (120, 120)])
        sage: p.xmin()
        100.0
        sage: p.ymin()
        100.0

    """
    if not isinstance(data[0], (list, tuple)):
        data = zip(range(len(data)),data)
    if plotjoined:
        P = line(data, **kwargs)
    else:
        P = point(data, **kwargs)
    return P

def networkx_plot(graph, pos=None, vertex_labels=True, vertex_size=300, vertex_colors=None,
                  edge_colors=None, graph_border=False, scaling_term=0.05, draw_edges=True):
    """
    Creates a graphics object ready to display a NetworkX graph.

    INPUT:
        graph -- a NetworkX graph
        pos -- an optional positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {   0: [-0.91679746, 0.88169588,],
                1: [ 0.47294849, 1.125     ,],
                2: [ 1.125     ,-0.12867615,],
                3: [ 0.12743933,-1.125     ,],
                4: [-1.125     ,-0.50118505,]   }
        vertex_labels -- determines whether labels for nodes are plotted
        vertex_size -- node size
        vertex_colors -- a dictionary specifying node colors: each key is a color recognized by
                        matplotlib, and each entry is a list of vertices.
        edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                        matplotlib, and each entry is a list of edges.
        scaling_term -- default is 0.05. if nodes are getting chopped off, increase; if graph
                        is too small, decrease. should be positive, but values much bigger than
                        1/8 won't be useful unless the nodes are huge

    EXAMPLES:
        sage: import networkx
        sage: D = networkx.dodecahedral_graph()
        sage: networkx_plot(D)

        sage: import networkx
        sage: from math import sin, cos, pi
        sage: P = networkx.petersen_graph()
        sage: d = {'#FF0000':[0,5], '#FF9900':[1,6], '#FFFF00':[2,7], '#00FF00':[3,8], '#0000FF':[4,9]}
        sage: pos_dict = {}
        sage: for i in range(5):
        ...    x = float(cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: for i in range(10)[5:]:
        ...    x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
        ...    y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
        ...    pos_dict[i] = [x,y]
        ...
        sage: networkx_plot(graph=P, vertex_colors=d, pos=pos_dict)

        sage: C = graphs.CubeGraph(5)
        sage: from sage.plot.plot import rainbow
        sage: R = rainbow(5)
        sage: edge_colors = {}
        sage: for i in range(5):
        ...    edge_colors[R[i]] = []
        sage: for u,v,l in C.edges():
        ...    for i in range(5):
        ...        if u[i] != v[i]:
        ...            edge_colors[R[i]].append((u,v,l))
        sage: networkx_plot(C.networkx_graph(), pos=C.get_pos(), edge_colors=edge_colors, vertex_labels=False, vertex_size=0)

    """
    g = Graphics()
    NGP = GraphicPrimitive_NetworkXGraph(graph, pos=pos, vertex_labels=vertex_labels, \
      vertex_size=vertex_size, vertex_colors=vertex_colors, edge_colors=edge_colors, \
      scaling_term=scaling_term, draw_edges=draw_edges)
    g._Graphics__objects.append(NGP)
    xmin = NGP._xmin
    xmax = NGP._xmax
    ymin = NGP._ymin
    ymax = NGP._ymax
    g.axes_range(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    if graph_border:
        from sage.plot.plot import line
        dx = (xmax - xmin)/10
        dy = (ymax - ymin)/10
        border = (line([( xmin - dx, ymin - dy), ( xmin - dx, ymax + dy ), ( xmax + dx, ymax + dy ), ( xmax + dx, ymin - dy ), ( xmin - dx, ymin - dy )], thickness=1.3))
        border.axes_range(xmin = (xmin - dx), xmax = (xmax + dx), ymin = (ymin - dy), ymax = (ymax + dy))
        g = g + border
    g.axes(False)
    return g

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

def to_mpl_color(c):
    """
    Convert a tuple or string to a matplotlib rgb color tuple.

    INPUT:
        c -- string or 3-tuple

    OUTPUT:
        3-tuple of floats between 0 and 1.

    EXAMPLES:
        sage: from sage.plot.plot import to_mpl_color
        sage: to_mpl_color('#fa0')
        (1.0, 0.66666666666666663, 0.0)
        sage: to_mpl_color('#ffffe1')
        (1.0, 1.0, 0.88235294117647056)
        sage: to_mpl_color('blue')
        (0.0, 0.0, 1.0)
        sage: to_mpl_color([1,1/2,1/3])
        (1.0, 0.5, 0.33333333333333331)
        sage: to_mpl_color([1,2,255])   # WARNING -- numbers are reduced mod 1!!
        (1.0, 0.0, 0.0)
    """
    if isinstance(c, Color):
        c = c.rgb()

    if isinstance(c, str):
        if len(c) > 0 and c[0] == '#':
            # it is some sort of html like color, e.g, #00ffff or #ab0
            h = c[1:]
            if len(h) == 3:
                h = '%s%s%s%s%s%s'%(h[0],h[0], h[1],h[1], h[2],h[2])
            elif len(h) != 6:
                raise ValueError, "color hex string (= '%s') must have length 3 or 6"%h
            return tuple([eval('0x%s'%h[i:i+2])/float(255) for i in [0,2,4]])
        else:
            from texture import colors
            try:
                return colors[c]
            except KeyError:
                raise ValueError, "unknown color '%s'"%c

    elif isinstance(c, (list, tuple)):
        c = list(c)
        if len(c) != 3:
            raise ValueError, "color tuple must have 3 entries, one for each RGB channel"
        for i in range(len(c)):
            s = float(c[i])
            if s != 1:
                s = modf(s)[0]
                if s < 0:
                    s += 1
            c[i] = s

    else:
        raise TypeError, "c must be a list, tuple, or string"

    return tuple(c)

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
        self._render(filename, dpi=dpi, figsize=self._figsize, **args)
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
        [(0.125000000000000*pi, 0.38268343236508978), (0.187500000000000*pi, 0.55557023301960218), (0.250000000000000*pi, 0.707106781186547...), (0.312500000000000*pi, 0.831469612302545...), (0.375000000000000*pi, 0.92387953251128674), (0.437500000000000*pi, 0.98078528040323043), (0.500000000000000*pi, 1.0), (0.562500000000000*pi, 0.98078528040323043), (0.625000000000000*pi, 0.92387953251128674), (0.687500000000000*pi, 0.831469612302545...), (0.750000000000000*pi, 0.70710678118654757), (0.812500000000000*pi, 0.55557023301960218), (0.875000000000000*pi, 0.3826834323650898...)]

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
