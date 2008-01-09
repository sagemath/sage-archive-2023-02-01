r"""
2D Plotting

Sage provides both Mathematica-style and Matlab-style plotting.

MATLAB-LIKE PLOTTING:
SAGE provides 2D plotting with an interface that is an exact
clone of Matlab (namely matplotlib).  For example,

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

MATHEMATICA-LIKE PLOTTING:
SAGE provides 2D plotting functionality with an interface inspired by
the interface for plotting in Mathematica.  The underlying rendering
is done using the matplotlib Python library.

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
    \item plot   -- plot of a function or other SAGE object (e.g., elliptic curve).
    \item parametric_plot
    \item polar_plot
    \item list_plot
    \item bar_chart
    \item contour_plot
    \item plot_vector_field
    \item matrix_plot
    \item graphics_array
\end{itemize}

The following misc Graphics functions are included:
\begin{itemize}
    \item Graphics
    \item is_Graphics
    \item rgbcolor
    \item hue
\end{itemize}

Type ? after each primitive in \sage for help and examples.

EXAMPLES:
We construct a plot involving several graphics objects:

    sage: G = plot(cos, -5, 5, thickness=5, rgbcolor=(0.5,1,0.5))
    sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))
    sage: P   # show it

Next we construct the reflection of the above polygon about the
$y$-axis by iterating over the qlist of first-coordinates of the first
graphic element of $P$ (which is the actual Polygon; note that $P$ is
a Graphics object, which consists of a single polygon):

    sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
    sage: Q   # show it

We combine together different graphics objects using "+":

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

Here is a pretty graph:
    sage: g = Graphics()
    sage: for i in range(60):
    ...    p = polygon([(i*cos(i),i*sin(i)), (0,i), (i,0)],\
    ...                rgbcolor=hue(i/40+0.4), alpha=0.2)
    ...    g = g + p
    ...
    sage: g.show(dpi=200, axes=False)

Another graph:
    sage: P = plot(lambda x: sin(x)/x, -4,4, rgbcolor=(0,0,1)) + \
    ...    plot(lambda x: x*cos(x), -4,4, rgbcolor=(1,0,0)) + \
    ...    plot(lambda x: tan(x),-4,4,rgbcolor=(0,1,0))
    ...
    sage: P.show(ymin=-pi,ymax=pi)

PYX EXAMPLES:
These are some examples of plots similar to some of the plots in the
PyX (http://pyx.sourceforge.net) documentation:

Symbolline:
    sage: y(x) = x*sin(x**2)
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
    sage: def f(x): return (x-3)*(x-5)*(x-7)+40
    sage: P = line([(2,0),(2,f(2))], rgbcolor=(0,0,0))
    sage: P += line([(8,0),(8,f(8))], rgbcolor=(0,0,0))
    sage: P += polygon([(2,0),(2,f(2))] + [(x, f(x)) for x in [2,2.1,..,8]] + [(8,0),(2,0)],  rgbcolor=(0.8,0.8,0.8))
    sage: P += text("$\\int_{a}^b f(x) dx$", (5, 20), fontsize=16, rgbcolor=(0,0,0))
    sage: P += plot(f, 1, 8.5, thickness=3)
    sage: P    # show the result

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

TODO:
    [] ability to change all properties of a graph easily, e.g.,
       the rgbcolor.
"""

############################################################################
#  Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com> and William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

from sage.structure.sage_object import SageObject

## IMPORTANT: Do *not* import matplotlib at module scope.  It takes a
## surprisingliy long time to initialize itself.  It's better if it is
## imported in functions, so it only gets started if it is actually
## going to be used.

DEFAULT_FIGSIZE=[6, 5]
DEFAULT_DPI = 100
EMBEDDED_MODE = False
DOCTEST_MODE = False
SHOW_DEFAULT = True

def show_default(default=None):
    """
    Set the default for showing plots using any plot commands.
    If called with no arguments, returns the current default.

    If this is True (the default) then any plot object when displayed
    will be displayed as an actual plot instead of text, i.e., the
    show command is not needed.

    EXAMPLES:
    The default starts out as True:
        sage: show_default()
        True

    We set it to False.
        sage: show_default(False)

    We see that it is False.
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

import random #for plot adaptive refinement
import os #for viewing and writing images
from colorsys import hsv_to_rgb #for the hue function
from math import sin, cos, modf, pi #for hue and polar_plot
from sage.structure.sage_object import SageObject

import sage.misc.misc

############### WARNING ###
# Try not to import any matplotlib stuff here -- matplotlib is
# slow to import.  (I did benchmarking and found that by not
# importing here, and instead importing when needed below, that
# SAGE startup times are much improved.)  - William
###############

#SAGE 2D Graphics Axes class:
from axes import Axes

def is_Graphics(x):
    """
    Return True if x is a Graphics object.

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

    Here we make a graphic of embeded isoceles triangles,
    coloring each one with a different color as we go:

        sage: h=10; c=0.4; p=0.1;
        sage: G = Graphics()
        sage: for x in srange(1,h+1):
        ...        l = [[0,x*sqrt(3)],[-x/2,-x*sqrt(3)/2],[x/2,-x*sqrt(3)/2],[0,x*sqrt(3)]]
        ...        G+=line(l,rgbcolor=hue(c + p*(x/h)))
        sage: G.show(figsize=[5,5])

    """

    def __init__(self):
        self.__xmin = -1
        self.__xmax = 1
        self.__ymin = -1
        self.__ymax = 1
        self.__fontsize = 8
        self.__show_axes = True
        self.__axes_color = (0, 0, 0)
        self.__axes_label_color = (0, 0, 0)
        self.__tick_color = (0, 0, 0)
        self.__tick_label_color = (0, 0, 0)
        self.__axes_width = 0.8
        self.__objects = []

    def range(self, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Set the ranges of the x and y axes.
        """
        self.xmin(xmin)
        self.xmax(xmax)
        self.ymin(ymin)
        self.ymax(ymax)

    def fontsize(self, s=None):
        """
        Set the font size of axes labels and tick marks.

        If called with no input, return the current fontsize.
        """
        if s is None:
            try:
                return self.__fontsize
            except AttributeError:
                self.__fontsize = 6
                return self.__fontsize
        self.__fontsize = s

    def axes(self, show=None):
        """
        Set whether or not the x and y axes are shown by default.

        If called with no input, return the current axes setting.
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
        """
        if c is None:
            try:
                return self.__axes_color
            except AttributeError:
                self.__axes_color = (0, 0, 0)
                return self.__axes_color
        self.__axes_color = c

    def axes_label(self, l=None):
        """
        Set the axes labels.

        If called with no input, return the current axes_label setting.
        """
        if l is None:
            try:
                return self.__axes_label
            except AttributeError:
                self.__axes_label = None
                return self.__axes_label
        self.__axes_label = l


    def axes_label_color(self, c=None):
        """
        Set the axes label color.

        If called with no input, return the current axes_label_color setting.
        """
        if c is None:
            try:
                return self.__axes_label_color
            except AttributeError:
                self.__axes_label_color = (0, 0, 0)
                return self.__axes_label_color
        self.__axes_label_color = c


    def axes_width(self, w=None):
        """
        Set the axes width.

        If called with no input, return the current axes_width setting.
        """
        if w is None:
            try:
                return self.__axes_width
            except AttributeError:
                self.__axes_width = True
                return self.__axes_width
        self.__axes_width = float(w)

    def tick_color(self, c=None):
        """
        Set the axes ticks color.

        If called with no input, return the current tick_color setting.
        """
        if c is None:
            try:
                return self.__tick_color
            except AttributeError:
                self.__tick_color = (0, 0, 0)
                return self.__tick_color
        self.__tick_color = c

    def tick_label_color(self, c=None):
        """
        Set the axes tick labels color.

        If called with no input, return the current tick_label_color setting.
        """
        if c is None:
            try:
                return self.__tick_label_color
            except AttributeError:
                self.__tick_label_color = (0, 0, 0)
                return self.__tick_label_color
        self.__tick_label_color = c



    def xmax(self, new=None):
        """
        sage: G = Graphics(); print G
        Graphics object consisting of 0 graphics primitives
        sage: G.xmax()
        1
        """
        if new is None:
            return self.__xmax
        self.__xmax = new

    def xmin(self, new=None):
        """
        sage: G = Graphics(); print G
        Graphics object consisting of 0 graphics primitives
        sage: G.xmin()
        -1
        """
        if new is None:
            return self.__xmin
        self.__xmin = new

    def ymax(self, new=None):
        """
        sage: G = Graphics(); print G
        Graphics object consisting of 0 graphics primitives
        sage: G.ymax()
        1
        """
        if new is None:
            return self.__ymax
        self.__ymax = new

    def ymin(self, new=None):
        """
        sage: G = Graphics(); print G
        Graphics object consisting of 0 graphics primitives
        sage: G.ymin()
        -1
        """
        if new is None:
            return self.__ymin
        self.__ymin = new

    def _repr_(self):
        if SHOW_DEFAULT:
            self.show()
            return ''
        else:
            return self.__str__()

    def __str__(self):
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
        if isinstance(other, int) and other == 0:
            return self
        raise TypeError

    def __add__(self, other):
        """
        If you have any Graphics object G1, you can
        always add any other amount of Graphics objects G2,G3,...
        to form a new Graphics object:
        G4 = G1 + G2 + G3

        EXAMPLES:
            sage: g1 = plot(abs(sqrt(x^3  - 1)), 1, 5)
            sage: g2 = plot(-abs(sqrt(x^3  - 1)), 1, 5)
            sage: g1 + g2  # displays the plot
        """
        if isinstance(other, int) and other == 0:
            return self
        if not isinstance(other, Graphics):
            raise TypeError, "other (=%s) must be a Graphics objects"%other
        g = Graphics()
        g.__xmin = min(self.__xmin, other.__xmin)
        g.__xmax = max(self.__xmax, other.__xmax)
        g.__ymin = min(self.__ymin, other.__ymin)
        g.__ymax = max(self.__ymax, other.__ymax)
        g.__objects = self.__objects + other.__objects
        return g

    def _arrow(self, xmin, ymin, xmax, ymax, options):
        self.__objects.append(GraphicPrimitive_Arrow(xmin, ymin, xmax, ymax, options))
        self._extend_axes(xmin, xmax, ymin, ymax)

    def _bar_chart(self, ind, datalist, xrange, yrange, options):
        self.__objects.append(GraphicPrimitive_BarChart(ind, datalist, options))
        self._extend_axes(xrange[0], xrange[1], yrange[0], yrange[1])

    def _circle(self, x, y, r, options):
        self.__objects.append(GraphicPrimitive_Circle(x, y, r, options))
        self._extend_axes(x+r, x-r, y+r, y-r)

    def _contour_plot(self, xy_data_array, xrange, yrange, options):
        self.__xmin = xrange[0]
        self.__xmax = xrange[1]
        self.__ymin = yrange[0]
        self.__ymax = yrange[1]
        self.__objects.append(GraphicPrimitive_ContourPlot(xy_data_array, xrange, yrange, options))

    def _disk(self, point, r, angle, options):
        xmin = point[0] - 2*r
        xmax = point[0] + 2*r
        ymin = point[1] - 2*r
        ymax = point[1] + 2*r
        self.__objects.append(GraphicPrimitive_Disk(point, r, angle, options))
        self._extend_axes(xmin, xmax, ymin, ymax)

    def _line(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Line(xdata, ydata, options))
        try:
            self._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass

    def _matrix_plot(self, xy_data_array, xrange, yrange, options):
        self.__xmin = xrange[0]
        self.__xmax = xrange[1]
        self.__ymin = yrange[0]
        self.__ymax = yrange[1]
        self.__objects.append(GraphicPrimitive_MatrixPlot(xy_data_array, xrange, yrange, options))
        #self._extend_axes(xrange[0], xrange[1], yrange[0], yrange[1])

    def _plot_field(self, xpos_array, ypos_array, xvec_array, yvec_array, xrange, yrange, options):
        self.__xmin = xrange[0]
        self.__xmax = xrange[1]
        self.__ymin = yrange[0]
        self.__ymax = yrange[1]
        self.__objects.append(GraphicPrimitive_PlotField(xpos_array, ypos_array, xvec_array, yvec_array, options))

    def _point(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Point(xdata, ydata, options))
        try:
            self._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass

    def _polygon(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Polygon(xdata, ydata, options))
        try:
            self._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass

    def _text(self, string, point, options):
        self.__objects.append(GraphicPrimitive_Text(string, point, options))
        xpad = 0.2*abs(point[0])
        ypad = 0.2*abs(point[1])
        self._extend_axes(point[0] - xpad, point[0] + xpad, point[1] - ypad, point[1] + ypad)

    def _extend_x_axis(self, x):
        xmin = self.__xmin
        xmax = self.__xmax
        if xmin is None or x < xmin:
            self.__xmin = x
        elif xmax is None or x > xmax:
            self.__xmax = x

    def _extend_y_axis(self, y):
        ymin = self.__ymin
        ymax = self.__ymax
        if ymin is None or y < ymin:
            self.__ymin = y
        elif ymax is None or y > ymax:
            self.__ymax = y

    def _extend_axes(self, xmin, xmax, ymin, ymax):
        self._extend_x_axis(xmin)
        self._extend_x_axis(xmax)
        self._extend_y_axis(ymin)
        self._extend_y_axis(ymax)

    def plot(self, *args, **kwds):
        return self

    def show(self, xmin=None, xmax=None, ymin=None, ymax=None,
             figsize=DEFAULT_FIGSIZE, filename=None,
             dpi=DEFAULT_DPI, axes=None, axes_label=None,frame=False,
             fontsize=None,
             **args):
        """
        Show this graphics image with the default image viewer.

        OPTIONAL INPUT:
            filename -- (default: None) string
            dpi -- dots per inch
            figsize -- [width, height] (same for square aspect)
            axes -- (default: True)
            fontsize -- positive integer
            frame -- (default: False) draw a MATLAB-like frame around the image

        EXAMPLES:
            sage: c = circle((1,1), 1, rgbcolor=(1,0,0))
            sage: c.show(xmin=-1, xmax=3, ymin=-1, ymax=3)

        To correct the apect ratio of certain graphics, it is necessary
        to show with a 'figsize' of square dimensions.

            sage: c.show(figsize=[5,5], xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can either turn off the drawing of the axes:

            sage: show(plot(sin,-4,4), axes=False)

        Or you can turn on the drawing of a frame around the plots:

            sage: show(plot(sin,-4,4), frame=True)

        """
        if DOCTEST_MODE:
            self.save(sage.misc.misc.SAGE_TMP + '/test.png',
                      xmin, xmax, ymin, ymax, figsize,
                    dpi=dpi, axes=axes, axes_label=axes_label,frame=frame)
            return
        if EMBEDDED_MODE:
            self.save(filename, xmin, xmax, ymin, ymax, figsize,
                    dpi=dpi, axes=axes, axes_label=axes_label,frame=frame)
            return
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self.save(filename, xmin, xmax, ymin, ymax, figsize, dpi=dpi, axes=axes,frame=frame, fontsize=fontsize)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(sage.misc.viewer.browser(), filename))

    def _prepare_axes(self, xmin, xmax, ymin, ymax):
        if self.__xmin is None:
            self.__xmin, self.__xmax, self.__ymin, self.__ymax = 0, 0, 0, 0

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
             dpi=DEFAULT_DPI, axes=None, axes_label=None, fontsize=None,
             frame=False, verify=True):
        """
        Save the graphics to an image file of type: PNG, PS, EPS, SVG, SOBJ,
        depending on the file extension you give the filename.
            Extension types can be: '.png', '.ps', '.eps', '.svg',
            and '.sobj' (for a SAGE object you can load later).

        EXAMPLES:
            sage: c = circle((1,1),1,rgbcolor=(1,0,0))
            sage: c.show(xmin=-1,xmax=3,ymin=-1,ymax=3)

            To correct the apect ratio of certain graphics, it is necessary
            to show with a 'figsize' of square dimensions.

            sage: c.show(figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)
        """
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

        self.axes_label(l=axes_label)
        #construct an Axes instance, see 'axes.py' for relevant code
        sage_axes = Axes(color=self.__axes_color, fontsize=self.__fontsize, axes_label=self.__axes_label,
                         axes_label_color=self.__axes_label_color, tick_color=self.__tick_color,
                         tick_label_color=self.__tick_label_color, linewidth=self.__axes_width)

        #adjust the xy limits and draw the axes:
        if not (contour or plotfield or matrixplot): #the plot is a 'regular' plot
            if frame: #add the frame axes
                xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
                axmin, axmax = xmin - 0.04*abs(xmax - xmin), xmax + 0.04*abs(xmax - xmin)
                aymin, aymax = ymin - 0.04*abs(ymax - ymin), ymax + 0.04*abs(ymax - ymin)
                subplot.set_xlim([axmin, axmax])
                subplot.set_ylim([aymin, aymax])
                #add a frame to the plot and possibly 'axes_with_no_ticks'
                sage_axes.add_xy_frame_axes(subplot, xmin, xmax, ymin, ymax,
                                        axes_with_no_ticks=axes, axes_label=axes_label)
            elif not frame and axes: #regular plot with regular axes
                xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
                subplot.set_xlim(xmin, xmax)
                subplot.set_ylim(ymin, ymax)
                sage_axes.add_xy_axes(subplot, xmin, xmax, ymin, ymax, axes_label=axes_label)
            else: #regular plot with no axes
                xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
                subplot.set_xlim(xmin, xmax)
                subplot.set_ylim(ymin, ymax)
        elif (contour or plotfield): #contour or field plot in self.__objects, so adjust axes accordingly
            xmin, xmax = self.__xmin, self.__xmax
            ymin, ymax = self.__ymin, self.__ymax
            subplot.set_xlim([xmin - 0.05*abs(xmax - xmin), xmax + 0.05*abs(xmax - xmin)])
            subplot.set_ylim([ymin - 0.05*abs(ymax - ymin), ymax + 0.05*abs(ymax - ymin)])
            if axes: #axes=True unless user specifies axes=False
                sage_axes.add_xy_frame_axes(subplot, xmin, xmax, ymin, ymax, axes_label=axes_label)
        else: #we have a 'matrix_plot' in self.__objects, so adjust axes accordingly
            xmin, xmax = self.__xmin, self.__xmax
            ymin, ymax = self.__ymin, self.__ymax
            subplot.set_xlim([xmin - 0.05*abs(xmax - xmin), xmax + 0.05*abs(xmax - xmin)])
            subplot.set_ylim([ymin - 0.05*abs(ymax - ymin), ymax + 0.05*abs(ymax - ymin)])
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
    def __init__(self, options):
        self.__options = options

    def _allowed_options(self):
        return {}


    def options(self):
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
        return "Graphics primitive"


class GraphicPrimitive_Arrow(GraphicPrimitive):
    """
    Primitive class that initializes the
    arrow graphics type
    """
    def __init__(self, xmin, ymin, xmax, ymax, options):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        return {'width':'How wide the entire arrow is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

    def _repr_(self):
        return "Arrow from (%s,%s) to (%s,%s) "%(self.xmin, self.ymin, self.xmax, self.ymax)

    def _render_on_subplot(self, subplot):
        options = self.options()
        width = float(options['width'])
        import matplotlib.patches as patches
        p = patches.FancyArrow(float(self.xmin), float(self.ymin), float(self.xmax), float(self.ymax),
                         width=width, length_includes_head=True)
        c = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(c)
        p.set_facecolor(c)
        subplot.add_patch(p)

#TODO: make bar_chart more general
class GraphicPrimitive_BarChart(GraphicPrimitive):
    """
    Primitive class that initializes the bar chart graphics primitive.
    """
    def __init__(self, ind, datalist, options):
        self.datalist = datalist
        self.ind = ind
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
        return {'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'width':'The width of the bars'}

    def _repr_(self):
        return "BarChart defined by a %s datalist "%(len(self.datalist))

    def _render_on_subplot(self, subplot):
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
    """
    def __init__(self, xdata, ydata, options):
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive.__init__(self, options)

    def _allowed_options(self):
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

    def _repr_(self):
        return "Line defined by %s points"%len(self)

    def __getitem__(self, i):
        return self.xdata[int(i)], self.ydata[int(i)]

    def __setitem__(self, i, point):
        i = int(i)
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        return len(self.xdata)

    def _render_on_subplot(self, subplot):
        import matplotlib.patches as patches
        options = dict(self.options())
        del options['alpha']
        del options['thickness']
        del options['rgbcolor']
        p = patches.lines.Line2D(self.xdata, self.ydata, **options)
        options = self.options()
        a = float(options['alpha'])
        p.set_alpha(a)
        p.set_linewidth(float(options['thickness']))
        p.set_color(to_mpl_color(options['rgbcolor']))
        subplot.add_line(p)

class GraphicPrimitive_Circle(GraphicPrimitive):
    """
    Primitive class that initializes the
    circle graphics type
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
                'contours':'Number of contour levels.'}

    def _repr_(self):
        return "ContourPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
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
            print "The possible color maps include: %s"%possibilities
            raise RuntimeError, "Color map %s not known"%cmap

        x0,x1 = float(self.xrange[0]), float(self.xrange[1])
        y0,y1 = float(self.yrange[0]), float(self.yrange[1])
        if fill:
            if contours is None:
                subplot.contourf(self.xy_data_array, cmap=cmap, extent=(x0,x1,y0,y1))
            else:
                subplot.contourf(self.xy_data_array, int(contours), cmap=cmap, extent=(x0,x1,y0,y1))
        else:
            if contours is None:
                subplot.contour(self.xy_data_array, cmap=cmap, extent=(x0,x1,y0,y1))
            else:
                subplot.contour(self.xy_data_array, int(contours), cmap=cmap, extent=(x0,x1,y0,y1))


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
            print "The possible color maps include: %s"%possibilities
            raise RuntimeError, "Color map %s not known"%cmap

        subplot.imshow(self.xy_data_array, cmap=cmap, interpolation='nearest', extent=(0,self.xrange[1],0,self.yrange[1]))

# Below is the base class that is used to make 'field plots'.
# Its implementation is motivated by Mathematica's 'PlotField'.
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
                'cmap':"""The colormap, one of (autumn, bone, cool, copper,
                gray, hot, hsv, jet, pink, prism, spring, summer, winter)"""}

    def _repr_(self):
        return "PlotField defined by a %s x %s vector grid"%(len(self.xpos_array), len(self.ypos_array))

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
            print "The possible color maps include: %s"%possibilities
            raise RuntimeError, "Color map %s not known"%cmap
        subplot.quiver(self.xpos_array, self.ypos_array, self.xvec_array, self.yvec_array)

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

    def _repr_(self):
        return "Point set defined by %s point(s)"%len(self.xdata)

    def __getitem__(self, i):
        return self.xdata[i], self.ydata[i]

    def _render_on_subplot(self,subplot):
        options = self.options()
        c = to_mpl_color(options['rgbcolor'])
        a = float(options['alpha'])
        s = int(options['pointsize'])
        faceted = options['faceted'] #faceted=True colors the edge of point
        subplot.scatter(self.xdata, self.ydata, s=s, c=c, alpha=a, faceted=faceted)

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
                'vertical_alignment': 'how to align vertically: top, center, bottom',
                'horizontal_alignment':'how to align horizontally: left, center, right'}

    def _render_on_subplot(self, subplot):
        options = self.options()
        c = options['rgbcolor']
        f = int(options['fontsize'])
        va = options['vertical_alignment']
        ha = options['horizontal_alignment']
        subplot.text(float(self.x), float(self.y), self.string, color=c, fontsize=f,
                        verticalalignment=va,horizontalalignment=ha)

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
        sage: pos = C.__get_pos__()
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
        sage: G.range(xmin=-1.1, xmax=2.2, ymin=0, ymax=3.25)
        sage: G.axes(False)
        sage: G.show()

    We color the edges and vertices of a Dodecahedral graph:
        sage: g = graphs.DodecahedralGraph()
        sage: g.show(edge_colors={(1.0, 0.8571428571428571, 0.0): g.edges()})

    """
    def __init__(self, graph, pos=None, vertex_labels=True, vertex_size=300, vertex_colors=None, edge_colors=None, scaling_term=0.05):
        self.__nxg = graph
        self.__vertex_size = vertex_size
        self.__vertex_labels = vertex_labels
        self.__vertex_colors = vertex_colors
        self.__edge_colors = edge_colors
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
                from random import random
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

######################################################################
#                                                                    #
#    Graphics Primitives Factories -- construct GraphicPrimitives    #
#                                                                    #
######################################################################
#
# The current method of writing a new Graphic Primitive
# involves writing a specific Factory for a given
# primitive, for example, 'GraphicPrimitive_circle' for 'circle'.
# This class should inherit from GraphicPrimitiveFactory,
# which should define any general Graphic Primitive attributes.
#
# As of now, the Graphic Primitive Factories, have
# only a __call__ method that deals with setting kwargs
# and coercing data into a correct form to present to
# one of the matplotlib functions.
#

class GraphicPrimitiveFactory:
    def __init__(self):
        # options for this specific graphics primitive.
        self.reset()

    def reset(self):
        # First the default options for all graphics primitives
        self.options = {'alpha':1,'thickness':1,'rgbcolor':(0,0,1)}
        self._reset()

    def _coerce(self, xdata, ydata):
        return to_float_list(xdata), to_float_list(ydata)

    def _graphic3d(self, *args, **kwds):
        """
        Return 3d version of this graphics primitive.

        We call this if the user tries to create a graphic but gives
        points (etc) in 3-space instead of in the plane.
        """
        raise NotImplementedError, "3d plotting of this primitive not yet implemented"

class GraphicPrimitiveFactory_arrow(GraphicPrimitiveFactory):
    def __call__(self, minpoint, maxpoint, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        xmin = float(minpoint[0])
        ymin = float(minpoint[1])
        xmax = float(maxpoint[0]) - xmin
        ymax = float(maxpoint[1]) - ymin
        return self._from_xdata_ydata(xmin, ymin, xmax, ymax, options=options)

#XXX: BarChart is a work in progress, only supports one datalist now.
class GraphicPrimitiveFactory_bar_chart(GraphicPrimitiveFactory):
    def __call__(self, datalist, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
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
        return self._from_xdata_ydata(ind, datalist, xrange, yrange, options=options)

class GraphicPrimitiveFactory_circle(GraphicPrimitiveFactory):
    def __call__(self, point, radius, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata((float(point[0]), float(point[1])),
                                       float(radius), options=options)

class GraphicPrimitiveFactory_contour_plot(GraphicPrimitiveFactory):
    def __call__(self, f, xrange, yrange, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        #should the xy_data_array be made right
        #when contour_plot is called?  Here we go:
        plot_points = int(options['plot_points'])
        xstep = abs(float(xrange[0]) - float(xrange[1]))/plot_points
        ystep = abs(float(yrange[0]) - float(yrange[1]))/plot_points
        xy_data_array = [[float(f(x, y)) for x in \
                          sage.misc.misc.xsrange(xrange[0], xrange[1], xstep)]
                         for y in sage.misc.misc.xsrange(yrange[0], yrange[1], ystep)]
        return self._from_xdata_ydata(xy_data_array, xrange, yrange, options=options)

class GraphicPrimitiveFactory_matrix_plot(GraphicPrimitiveFactory):
    def __call__(self, mat, **kwds):
        from sage.matrix.all import is_Matrix
        from matplotlib.numerix import array
        if not is_Matrix(mat) or (isinstance(mat, (list, tuple)) and isinstance(mat[0], (list, tuple))):
            raise TypeError, "mat must be of type Matrix or a two dimensional array"
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        if is_Matrix(mat):
            xrange = (0, mat.ncols())
            yrange = (0, mat.nrows())
        else:
            xrange = (0, len(mat[0]))
            yrange = (0, len(mat))
        xy_data_array = [array(r, dtype=float) for r in mat]
        return self._from_xdata_ydata(xy_data_array, xrange, yrange, options=options)


class GraphicPrimitiveFactory_plot_field(GraphicPrimitiveFactory):
    def __call__(self, (f, g), xrange, yrange, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        #big list loop, again, should this be done here?:
        plot_points = int(options['plot_points'])
        xstep = abs(float(xrange[0]) - float(xrange[1]))/plot_points
        ystep = abs(float(yrange[0]) - float(yrange[1]))/plot_points
        Lpx,Lpy,Lcx,Lcy = [],[],[],[]
        for x in sage.misc.misc.xsrange(xrange[0], xrange[1], xstep):
            for y in sage.misc.misc.xsrange(yrange[0], yrange[1], ystep):
                Lpx.append(float(x))
                Lpy.append(float(y))
                Lcx.append(float(f(x)))
                Lcy.append(float(g(y)))
        return self._from_xdata_ydata(Lpx, Lpy, Lcx, Lcy, xrange, yrange, options=options)

class GraphicPrimitiveFactory_disk(GraphicPrimitiveFactory):
    def __call__(self, point, radius, angle, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata((float(point[0]), float(point[1])),float(radius),
                    (float(angle[0]), float(angle[1])), options=options)

class GraphicPrimitiveFactory_text(GraphicPrimitiveFactory):
    def __call__(self, string, point, **kwds):
        if len(point) == 3:
            from sage.plot.plot3d.shapes2 import text3d
            return text3d(string, point, **kwds)
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata(string, (float(point[0]), float(point[1])), options=options)

class GraphicPrimitiveFactory_points(GraphicPrimitiveFactory):
    def __call__(self, xdata, ydata, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata(xdata, ydata, options=options)

# WARNING: The below GraphicPrimitiveFactory_from_point_list
# class can potentially be very slow for large point sets.
#
# It exists because it provides the following functionality:
# Allows user to give as input to the function 'point'
# a list of (x,y) values at which to plot points, coloring
# each one a different color if needed.  From this input list we then
# loop through it, first coercing all the values the floats
# and then forming two new list that consist of all then
# x-values in one list and all the y-values in another list.
# This is needed to be done because that is how the input is
# taken in the matplotlib function 'scatter'.

class GraphicPrimitiveFactory_from_point_list(GraphicPrimitiveFactory):
    def __call__(self, points, coerce=True, **kwds):
        try:
            return points.plot(**kwds)
        except AttributeError:
            pass
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v

        if not isinstance(points, (list,tuple)) or \
           (isinstance(points,(list,tuple)) and len(points) <= 3 and not
            isinstance(points[0], (list,tuple))):
            try:
                points = [[float(z) for z in points]]
            except TypeError:
                pass

        try:
            if len(points) > 0 and len(points[0]) == 3:
                return self._graphic3d()(points, coerce=coerce, **kwds)
        except (AttributeError, TypeError):
            pass
        xdata = []
        ydata = []
        if coerce:
            xdata = [float(z[0]) for z in points]
            ydata = [float(z[1]) for z in points]
        else:
            xdata = [z[0] for z in points]
            ydata = [z[1] for z in points]

        return self._from_xdata_ydata(xdata, ydata, True, options=options)


class ArrowFactory(GraphicPrimitiveFactory_arrow):
    """

    An arrow from (xmin, ymin) to (xmax, ymax).

    EXAMPLES:

    A straight, black arrow
       sage: arrow((1, 1), (3, 3))

    Make a red arrow:
       sage: arrow((-1, -1), (2, 3), rgbcolor=(1,0,0))

    You can change the width of an arrow:
        sage: arrow((1, 1), (3, 3), width=0.05)
    """
    def _reset(self):
        self.options={'width':0.02,'rgbcolor':(0, 0, 1)}

    def _repr_(self):
        return "type arrow? for help and examples"

    def _from_xdata_ydata(self, xmin, ymin, xmax, ymax, options):
        g = Graphics()
        g._arrow(xmin, ymin, xmax, ymax, options=options)
        return g

#an unique arrow instance
arrow = ArrowFactory()

class BarChartFactory(GraphicPrimitiveFactory_bar_chart):
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
    def _reset(self):
        self.options={'width':0.5,'rgbcolor':(0, 0, 1)}

    def _repr_(self):
        return "type bar_chart? for help and examples"

    def _from_xdata_ydata(self, ind, datalist, xrange, yrange, options):
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

#an unique bar_chart instance
bar_chart = BarChartFactory()


class CircleFactory(GraphicPrimitiveFactory_circle):
    """
    Return a circle at a point = (x,y) with radius = r.
    Type circle.options to see all options

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
    to show with a 'figsize' of square dimensions.

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
    """
    def _reset(self):
        self.options={'alpha':1,'fill':False,'thickness':1,'rgbcolor':(0, 0, 1)}

    def _repr_(self):
        return "type circle? for help and examples"

    def _from_xdata_ydata(self, point, r, options):
        g = Graphics()
        g._circle(float(point[0]), float(point[1]), float(r), options)
        return g


#an unique circle instance
circle = CircleFactory()

class ContourPlotFactory(GraphicPrimitiveFactory_contour_plot):
    r"""

    \code{contour_plot} takes a function of two variables, f(x,y)
    and plots contour lines of the function over the specified
    xrange and yrange as demonstrated below.

      contour_plot(f, (xmin, xmax), (ymin, ymax), ...)

    INPUT:
        f -- a function of two variables
        (xmin, xmax) -- 2-tuple, the range of x values
        (ymin, ymax) -- 2-tuple, the range of y values
    The following inputs must all be passed in as named parameters:
        plot_points  -- integer (default: 25); number of points to plot
                        in each direction of the grid
        fill         -- bool (default: True), whether to color in the area
                        between contour lines
        cmap         -- string (default: 'gray'), the color map to use:
                        autumn, bone, cool, copper, gray, hot, hsv,
                        jet, pink, prism, spring, summer, winter
        contours     -- integer (default: None), the number of contour
                        lines to draw.  If None, determined automatically,
                        and usually about 5.


    EXAMPLES:

    Here we plot a simple funtion of two variables:
        sage: def f(x,y):
        ...       return cos(x^2 + y^2)
        sage: contour_plot(f, (-4, 4), (-4, 4))


    Here we change the ranges and add some options:
        sage: def h(x,y):
        ...       return (x^2)*cos(x*y)
        sage: contour_plot(h, (-10, 5), (-5, 5), fill=False, plot_points=100)


    An even more complicated plot.
        sage: def f(x,y):
        ...       return sin(x^2 + y^2)*cos(x)*sin(y)
        sage: contour_plot(f, (-4, 4), (-4, 4),plot_points=100)
    """
    def _reset(self):
        self.options={'plot_points':25, 'fill':True, 'cmap':'gray', 'contours':None}

    def _repr_(self):
        return "type contour_plot? for help and examples"

    def _from_xdata_ydata(self, xy_data_array, xrange, yrange, options):
        g = Graphics()
        g._contour_plot(xy_data_array, xrange, yrange, options)
        return g

#unique contour_plot instance
contour_plot = ContourPlotFactory()

class LineFactory(GraphicPrimitiveFactory_from_point_list):
    r"""
    Create the line through the given list of points.

    Type line.options for a dictionary of the default options for
    lines.  You can change this to change the defaults for all future
    lines.  Use line.reset() to reset to the default options.

    INPUT:
    \begin{verbatim}
        alpha -- How transparent the line is
        thickness -- How thick the line is
        rgbcolor -- The color as an rgb tuple
        hue -- The color given as a hue
        Any MATLAB/MATPLOTLIB line option may also be passed in.  E.g.,
        linestyle -- The style of the line, which is one of
                  '--' (dashed), '-.' (dash dot), '-' (solid),
                  'steps', ':' (dotted)
        marker  -- "'0' (tickleft), '1' (tickright), '2' (tickup), '3' (tickdown),
                   '' (nothing), ' ' (nothing), '+' (plus), ',' (pixel), '.' (point),
                   '1' (tri_down), '3' (tri_left), '2' (tri_up), '4' (tri_right),
                   '<' (triangle_left), '>' (triangle_right), 'None' (nothing),
                   'D' (diamond), 'H' (hexagon2), '_' (hline), '^' (triangle_up),
                   'd' (thin_diamond), 'h' (hexagon1), 'o' (circle), 'p' (pentagon),
                   's' (square), 'v' (triangle_down), 'x' (x), '|' (vline)"
       markersize -- the size of the marker in points
       markeredgecolor -- the markerfacecolor can be any color arg
       markeredgewidth -- the size of the markter edge in points
    \end{verbatim}


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
        sage: vals = E.Lseries().values_along_line(1-I, 1+10*I, 100) # critical line
        sage: L = [(z[1].real(), z[1].imag()) for z in vals]
        sage: line(L, rgbcolor=(3/4,1/2,5/8))

    A red, blue, and green "cool cat":

        sage: G = plot(-cos(x), -2, 2, thickness=5, rgbcolor=(0.5,1,0.5))
        sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))
        sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
        sage: G + P + Q   # show the plot
    """
    def _reset(self):
        self.options = {'alpha':1,'rgbcolor':(0,0,1),'thickness':1}

    def _repr_(self):
        return "type line? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g._Graphics__objects.append(GraphicPrimitive_Line(xdata, ydata, options))
        try:
            g._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass
        return g

    def _graphic3d(self):
        from sage.plot.plot3d.shapes2 import line3d
        return line3d

# unique line instance
line = LineFactory()

class MatrixPlotFactory(GraphicPrimitiveFactory_matrix_plot):
    r"""
    A plot of a given matrix or 2D array.

    Each (ith, jth) matrix element is given a different
    color value depending on its relative size compared
    to the other elements in the matrix.

    The tick marks drawn on the frame axes denote the
    (ith, jth) element of the matrix.

    EXAMPLES:

    A matrix over ZZ colored with different grey levels:

        sage: matrix_plot(matrix([[1,3,5,1],[2,4,5,6],[1,3,5,7]]))

    Here we make a random matrix over RR and use cmap='hsv'
    to color the matrix elements different RGB colors:

        sage: matrix_plot(random_matrix(RDF, 50), cmap='hsv')

    Another random plot, but over GF(389):
        sage: matrix_plot(random_matrix(GF(389), 10), cmap='Oranges')
    """
    def _reset(self):
        self.options={'cmap':'gray'}

    def _repr_(self):
        return "type matrix_plot? for help and examples"

    def _from_xdata_ydata(self, xy_data_array, xrange, yrange, options):
        g = Graphics()
        g._matrix_plot(xy_data_array, xrange, yrange, options)
        return g

#unique matrix_plot instance
matrix_plot = MatrixPlotFactory()


# Below is the base class that is used to make 'plot_vector_field'.
# Its implementation is motivated by Mathematica's 'PlotVectorField'.
# TODO: make class similiar to this one to implement:
# 'plot_gradient_field' and 'plot_hamiltonian_field'
class PlotFieldFactory(GraphicPrimitiveFactory_plot_field):
    r"""

    \code{plot_field} takes two functions of one variable, (f(x), g(y))
    and plots vector arrows of the function over the specified
    xrange and yrange as demonstrated below.
    plot_field((f, g), (xmin, xmax), (ymin, ymax))

    EXAMPLES:

    Plot the vector field of sin and cos:
    sage: vf1 = plot_vector_field((lambda x:sin(x), lambda y:cos(y)), (-3,3), (-3,3))

    """
    def _reset(self):
        self.options={'plot_points':20, 'cmap':'gray'}

    def _repr_(self):
        return "type contour_plot? for help and examples"

    def _from_xdata_ydata(self, xpos_array, ypos_array, xvec_array, yvec_array, xrange, yrange, options):
        g = Graphics()
        g._plot_field(xpos_array, ypos_array, xvec_array, yvec_array, xrange, yrange, options)
        return g

#unique plot_vector_field instance
plot_vector_field = PlotFieldFactory()


class DiskFactory(GraphicPrimitiveFactory_disk):
    """

    A disk at a point = (x,y) with radius = r
    spanning (in radians) angle=(rad1, rad2)
    Type disk.options to see all options

    EXAMPLES:
    Make some dangerous disks:

        sage: bl = disk((0.0,0.0), 1, (pi, 3*pi/2), rgbcolor=(1,1,0))
        sage: tr = disk((0.0,0.0), 1, (0, pi/2), rgbcolor=(1,1,0))
        sage: tl = disk((0.0,0.0), 1, (pi/2, pi), rgbcolor=(0,0,0))
        sage: br = disk((0.0,0.0), 1, (3*pi/2, 2*pi), rgbcolor=(0,0,0))
        sage: P  = tl+tr+bl+br
        sage: P.show(figsize=(4,4),xmin=-2,xmax=2,ymin=-2,ymax=2)

    """
    def _reset(self):
        self.options={'alpha':1,'fill':True,'rgbcolor':(0,0,1),'thickness':0}

    def _repr_(self):
        return "type disk? for help and examples"

    def _from_xdata_ydata(self, point, r, angle, options):
        g = Graphics()
        g._disk(point, r, angle, options)
        return g

#an unique disk instance
disk = DiskFactory()

class PointFactory(GraphicPrimitiveFactory_from_point_list):
    """

    A point of size 'pointsize' defined by point = (x,y)
    Type point.options to see all options. point takes either
    a single tuple of coordinates or a list of tuples.

    EXAMPLES:
        A purple point from a single tuple or coordinates:
        sage: point((0.5, 0.5), rgbcolor=hue(0.75))

        Here are some random larger red points, given as a list of tuples
        sage: point(((0.5, 0.5), (1, 2), (0.5, 0.9), (-1, -1)), rgbcolor=hue(1), pointsize=30)

    """
    def _reset(self):
        self.options = {'alpha':1,'pointsize':10,'faceted':False,'rgbcolor':(0,0,1)}

    def _repr_(self):
        return "type point? for options help"

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g._Graphics__objects.append(GraphicPrimitive_Point(xdata, ydata, options))
        try:
            g._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass
        return g

# unique point instance
point = PointFactory()
points = point


class PolygonFactory(GraphicPrimitiveFactory_from_point_list):
    """
    Type polygon.options for a dictionary of the default
    options for polygons.  You can change this to change
    the defaults for all future polygons.  Use polygon.reset()
    to reset to the default options.

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

    AUTHORS:
        -- David Joyner (2006-04-14): the long list of examples above.

    """
    def _reset(self):
        self.options={'alpha':1,'rgbcolor':(0,0,1),'thickness':0}

    def _repr_(self):
        return "SAGE polygon; type polygon? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g._Graphics__objects.append(GraphicPrimitive_Polygon(xdata, ydata, options))
        try:
            g._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass
        return g

# unique polygon instance
polygon = PolygonFactory()

class PlotFactory(GraphicPrimitiveFactory):
    r"""
    Use plot by writing

        plot(X, ...)

    where X is a Sage object (or list of Sage objects) that either is
    callable and returns numbers that can be coerced to floats, or has
    a plot method that returns a GraphicPrimitive object.

    Type plot.options for a dictionary of the default
    options for plots.  You can change this to change
    the defaults for all future plots.  Use plot.reset()
    to reset to the default options.

    PLOT OPTIONS:
    The plot options are

        plot_points -- the number of points to initially plot before
                       doing adaptive refinement
        plot_division -- the maximum number points including those
                       computed during adaptive refinement
        max_bend      -- parameter that affects adaptive refinement

        xmin -- starting x value
        xmax -- ending x value

    APPEARANCE OPTIONS:
    The following options affect the appearance of the line through the points
    on the graph of X (these are the same as for the line function):

    INPUT:
    \begin{verbatim}
        alpha -- How transparent the line is
        thickness -- How thick the line is
        rgbcolor -- The color as an rgb tuple
        hue -- The color given as a hue
        Any MATLAB/MATPLOTLIB line option may also be passed in.  E.g.,
        linestyle -- The style of the line, which is one of
                  '--' (dashed), '-.' (dash dot), '-' (solid),
                  'steps', ':' (dotted)
        marker  -- "'0' (tickleft), '1' (tickright), '2' (tickup), '3' (tickdown),
                   '' (nothing), ' ' (nothing), '+' (plus), ',' (pixel), '.' (point),
                   '1' (tri_down), '3' (tri_left), '2' (tri_up), '4' (tri_right),
                   '<' (triangle_left), '>' (triangle_right), 'None' (nothing),
                   'D' (diamond), 'H' (hexagon2), '_' (hline), '^' (triangle_up),
                   'd' (thin_diamond), 'h' (hexagon1), 'o' (circle), 'p' (pentagon),
                   's' (square), 'v' (triangle_down), 'x' (x), '|' (vline)"
       markersize -- the size of the marker in points
       markeredgecolor -- the markerfacecolor can be any color arg
       markeredgewidth -- the size of the marker edge in points
    \end{verbatim}

    Note that this function does NOT simply sample equally spaced
    points between xmin and xmax.  Instead it computes equally spaced
    points and add small perturbations to them.  This reduces the
    possibility of, e.g., sampling sin only at multiples of $2\pi$,
    which would yield a very misleading graph.

    EXAMPLES:
    We plot the sin function:
        sage: P = plot(sin, 0,10); print P
        Graphics object consisting of 1 graphics primitive
        sage: len(P)     # number of graphics primitives
        1
        sage: len(P[0])  # how many points were computed
        200
        sage: P          # render

        sage: P = plot(sin, 0,10, plot_points=10); print P
        Graphics object consisting of 1 graphics primitive
        sage: len(P[0])  # random output
        80
        sage: P          # render

    We plot several functions together by passing a list
    of functions as input:
       sage: plot([sin(n*x) for n in [1..4]], 0, pi)


    The function $\sin(1/x)$ wiggles wildtly near $0$, so the
    first plot below won't look perfect.  Sage adapts to this
    and plots extra points near the origin.
       sage: plot(sin(1/x), -1, 1)

    The \code{plot_points} option, you can increase the number
    of sample points, to obtain a more accurate plot.
       sage: plot(sin(1/x), -1, 1, plot_points=1000)

    The actual sample points are slightly randomized, so the above
    plots may look slightly different each time you draw them.

    We draw the graph of an elliptic curve as the union
    of graphs of 2 functions.
        sage: def h1(x): return abs(sqrt(x^3  - 1))
        sage: def h2(x): return -abs(sqrt(x^3  - 1))
        sage: P = plot([h1, h2], 1,4)    # slightly random output because of random sampling
        Graphics object consisting of 2 graphics primitives
        sage: P          # show the result

    We can also directly plot the elliptic curve:
        sage: E = EllipticCurve([0,-1])
        sage: plot(E, 1, 4, rgbcolor=hue(0.6))

    We can change the line style to one of '--' (dashed), '-.' (dash dot),
    '-' (solid), 'steps', ':' (dotted):
        sage: plot(sin(x), 0, 10, linestyle='-.')
    """
    def _reset(self):
        o = self.options
        o['plot_points'] = 200
        o['plot_division'] = 1000
        o['max_bend'] = 0.1
        o['rgbcolor'] = (0,0,1)

    def _repr_(self):
        return "plot; type plot? for help and examples."

    def __call__(self, funcs, *args, **kwds):
        do_show = False
        if kwds.has_key('show') and kwds['show']:
            do_show = True
            del kwds['show']
        if hasattr(funcs, 'plot'):
            G = funcs.plot(*args, **kwds)
        # if we are using the generic plotting method
        else:
            n = len(args)
            # if there is one extra arg, then it had better be a tuple
            if n == 1:
                G = self._call(funcs, *args, **kwds)
            elif n == 2:
            # if ther eare two extra args, then pull them out and pass them as a tuple
                xmin = args[0]
                xmax = args[1]
                args = args[2:]
                G = self._call(funcs, (xmin, xmax), *args, **kwds)
            else:
                print "there were %s extra arguments (besides %s)" % (n, funcs)
        if do_show:
            G.show()
        return G

    def _call(self, funcs, xrange, parametric=False,
              polar=False, label='', **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v

        #parametric_plot will be a list or tuple of two functions (f,g)
        #and will plotted as (f(x), g(x)) for all x in the given range
        if parametric:
            if len(funcs) == 3:
                raise ValueError, "use parametric_plot3d for parametric plots in 3d dimensions."
            elif len(funcs) == 2:
                # 2d
                f,g = funcs
            else:
                raise ValueError, "parametric plots only implemented in 2 and 3 dimensions."

        #or we have only a single function to be plotted:
        else:
            f = funcs

        plot_points = int(options['plot_points'])
        del options['plot_points']
        x, data = var_and_list_of_values(xrange, plot_points)
        data = list(data)
        xmin = data[0]
        xmax = data[-1]

        #check to see if funcs is a list of functions that will
        #be all plotted together.
        if isinstance(funcs, (list, tuple)) and not parametric:
            G = Graphics()
            for i in range(0, len(funcs)):
                G += plot(funcs[i], (xmin, xmax), polar=polar, **kwds)
            return G

        delta = (xmax - xmin) / plot_points
        dd = delta

        exceptions = 0; msg=''
        for i in range(plot_points):
            xi = xmin + i*delta
            if i < plot_points:
                xi += delta*random.random()
                if xi > xmax:
                    xi = xmax
            else:
                xi = xmax  # guarantee that we get the last point.

            try:
                y = f(xi)
                data[i] = (float (xi), float(y))
            except (ZeroDivisionError, TypeError, ValueError), msg:
                sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, x),1)
                exceptions += 1

        # adaptive refinement
        i, j = 0, 0
        max_bend = float(options['max_bend'])
        del options['max_bend']
        plot_division = int(options['plot_division'])
        del options['plot_division']
        while i < len(data) - 1:
            if abs(data[i+1][1] - data[i][1]) > max_bend:
                x = (data[i+1][0] + data[i][0])/2
                try:
                    y = float(f(x))
                    data.insert(i+1, (x, y))
                except (ZeroDivisionError, TypeError, ValueError), msg:
                    sage.misc.misc.verbose("%s\nUnable to compute f(%s)"%(msg, x),1)
                    exceptions += 1

                j += 1
                if j > plot_division:
                    break
            else:
                i += 1

        if (len(data) == 0 and exceptions > 0) or exceptions > 10:
            print "WARNING: When plotting, failed to evaluate function at %s points."%exceptions
            print "Last error message: '%s'"%msg
        if parametric:
            data = [(fdata, g(x)) for x, fdata in data]
        if polar:
            data = [(y*cos(x), y*sin(x)) for x, y in data]
        G = line(data, coerce=False, **options)

        # Label?
        if label:
            label = '  '+str(label)
            G += text(label, data[-1], horizontal_alignment='left',
                      vertical_alignment='center')

        return G

# unique plot instance
plot = PlotFactory()


class TextFactory(GraphicPrimitiveFactory_text):
    """
    text(txt, point, **kwds):

    Returns a 2d or 3d text graphics object at the point (x,y)

    Type text.options for a dictionary of options for 2d text.  The 3d options
    are as for other 3d graphics objects (i.e., mainly just rgbcolor at present).

    2D OPTIONS:
        fontsize -- How big the text is
        rgbcolor -- The color as an rgb tuple
        hue -- The color given as a hue
        vertical_alignment -- how to align vertically: top, center, bottom
        horizontal_alignment -- how to align horizontally: left, center, right

    3D OPTIONS:
        rgbcolor -- the color of the text

    EXAMPLES:
    Some 2d text:
        sage: text("SAGE is really neat!!",(2,12))

    The same text, but in 3d:
        sage: text("SAGE is really neat!!",(2,12,1))

    The same text in larger font and colored red:
        sage: text("SAGE is really neat!!",(2,12),fontsize=20,rgbcolor=(1,0,0))

    And in 3d in two places:
        sage: text("SAGE is...",(2,12,1), rgbcolor=(1,0,0)) + text("quite powerful!!",(4,10,0), rgbcolor=(0,0,1))

    You can also align 2d text differently:
        sage: t1 = text("Hello",(1,1), vertical_alignment="top")
        sage: t2 = text("World", (1,0.5), horizontal_alignment="left")
        sage: t1 + t2   # render the sume
    """
    def _reset(self):
        self.options = {'fontsize':10, 'rgbcolor':(0,0,1),
                        'horizontal_alignment':'center',
                        'vertical_alignment':'center'}

    def _repr_(self):
        return "type text? for help and examples"

    def _from_xdata_ydata(self, string, point, options):
        g = Graphics()
        g._text(string, point, options)
        return g

# unique text instance
text = TextFactory()


########## misc functions ###################

def parametric_plot(funcs, tmin, tmax, **kwargs):
    """
    parametric_plot takes two functions as a list or a tuple and make
    a plot with the first function giving the x coordinates and the
    second function giving the y coordinates.

    INPUT:
        funcs -- 2 or 3-tuple of functions
        tmin -- start value of t
        tmax -- end value of t
        other options -- passed to plot.

    EXAMPLE:
    We draw a 2d parametric plot:
        sage: t = var('t')
        sage: parametric_plot( (sin(t), sin(2*t)), 0, 2*pi, rgbcolor=hue(0.6) )

    We draw a 3d parametric plot:
        sage: parametric_plot3d( (5*cos(x), 5*sin(x), x), (-12, 12), plot_points=150, color="red")
    """
    return plot(funcs, tmin, tmax, parametric=True, **kwargs)

def polar_plot(funcs, xmin, xmax, **kwargs):
    """
    polar_plot takes a single function or a list or tuple of functions
    and plots them parametrically in the given range.

    EXAMPLES:
    Here is a blue 8-leaved petal:
        sage: polar_plot(lambda x:sin(5*x)^2, 0, 2*pi, rgbcolor=hue(0.6))

    A red figure-8:
        sage: polar_plot(lambda x:abs(sqrt(1 - sin(x)^2)), 0, 2*pi, rgbcolor=hue(1.0))

    A green limacon of Pascal:
        sage: polar_plot(lambda x:2 + 2*cos(x), 0, 2*pi, rgbcolor=hue(0.3))

    """
    return plot(funcs, xmin, xmax, polar=True, **kwargs)

def list_plot(data, plotjoined=False, **kwargs):
    """
    list_plot takes a single list of data, in which case it forms a
    list of tuples (i,di) where i goes from 0 to len(data)-1 and di is
    the ith data value, and puts points at those tuple values.

    list_plot also takes a list of tuples (dxi, dyi) where dxi is the
    ith data representing the x-value, and dyi is the ith y-value if
    plotjoined=True, then a line spanning all the data is drawn
    instead.

    EXAMPLES:
        sage: list_plot([i^2 for i in range(5)])

    Here are a bunch of random red points:
        sage: r = [(random(),random()) for _ in range(20)]
        sage: list_plot(r,rgbcolor=(1,0,0))

    This gives all the random points joined in a purple line:
        sage: list_plot(r, plotjoined=True, rgbcolor=(1,0,1))
    """
    if not isinstance(data[0], (list, tuple)):
        data = zip(range(len(data)),data)
    if plotjoined:
        P = line(data, **kwargs)
    else:
        P = point(data, **kwargs)
    return P

def networkx_plot(graph, pos=None, vertex_labels=True, vertex_size=300, vertex_colors=None,
                  edge_colors=None, graph_border=False, scaling_term=0.05):
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
        sage: networkx_plot(C._nxg, pos=C.__get_pos__(), edge_colors=edge_colors, vertex_labels=False, vertex_size=0)
    """
    g = Graphics()
    NGP = GraphicPrimitive_NetworkXGraph(graph, pos=pos, vertex_labels=vertex_labels, vertex_size=vertex_size, vertex_colors=vertex_colors, edge_colors=edge_colors, scaling_term=scaling_term)
    g._Graphics__objects.append(NGP)
    xmin = NGP._xmin
    xmax = NGP._xmax
    ymin = NGP._ymin
    ymax = NGP._ymax
    g.range(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    if graph_border:
        from sage.plot.plot import line
        dx = (xmax - xmin)/10
        dy = (ymax - ymin)/10
        border = (line([( xmin - dx, ymin - dy), ( xmin - dx, ymax + dy ), ( xmax + dx, ymax + dy ), ( xmax + dx, ymin - dy ), ( xmin - dx, ymin - dy )], thickness=1.3))
        border.range(xmin = (xmin - dx), xmax = (xmax + dx), ymin = (ymin - dy), ymax = (ymax + dy))
        g = g + border
    g.axes(False)
    return g

def to_float_list(v):
    return [float(x) for x in v]

def to_mpl_color(c):
    c = list(c)
    for i in range(len(c)):
        s = float(c[i])
        if s != 1:
            s = modf(s)[0]
            if s < 0:
                s += 1
        c[i] = s
    return tuple(c)

def hue(h, s=1, v=1):
    """
      hue(h,s=1,v=1) where 'h' stands for hue,
      's' stands for saturation, 'v' stands for value.
      hue returns a list of rgb intensities (r, g, b)
      All values are in range 0 to 1.

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
    GraphicsArray takes a (m x n) list of lists of graphics
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
            frame -- (default: False) draw a MATLAB-like frame around the image
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
    Ten you can type either: \code{ga.show()} or \code{ga.save()}.

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
    cases where matplotlib is just plain schizophrenic- for an example, do

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
    Given an integer n, returns a list of colors, represented in HTML hex,
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
        plot_points -- integer >= 2 (the endpoints)
        v -- (v0, v1) or (var, v0, v1); if the former return
             the range of values between v0 and v1 taking
             plot_points steps; if var is given, also return var.

    OUTPUT:
        var -- a variable or None
        list -- a list of floats
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
        step = (b-a)/float(plot_points)
        values = [a + step*i for i in xrange(plot_points)]
        # want to make sure that we plot exactly as many points as requested
#         rng.append(b)
        return var, values

# def float_range(a, b, step):
#     """
#     Returns a list of floating point numbers from a to b with the
#     given step
#     """
#     (a,b,step) = (float(a),float(b),float(step))
#     v = [a]
#     w = a + step
#     while w < b:
#         v.append(w)
#         w += step
#     if w < b:
#         v.append(b)
#     return v
