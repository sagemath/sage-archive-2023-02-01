r"""
2D Plotting

SAGE provides 2-d plotting functionality with an interface inspired by
the interface for plotting in Mathematica.  The underlying rendering
is mostly implemented using the matplotlib Python library.

The following graphics primitives are supported:
\begin{itemize}
    \item line   -- a line determined by a sequence of points (this need not be straight!)
    \item circle -- a circle with given radius
    \item disk   -- a filled disk
    \item point  -- a point
    \item text   -- some text
    \item polygon -- a filled polygon
    \item plot   -- plot of a function or other SAGE object (e.g., elliptic curve).
    \item parametric_plot
    \item list_plot
\end{itemize}

Type ? after each primitive in \sage for help and examples.

EXAMPLES:
We construct a plot involving several graphics objects:

    sage: G = plot(cos, -5, 5, thickness=5, rgbcolor=(0.5,1,0.5))
    sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))

Next we construct the reflection of the above polygon about the
$y$-axis by iterating over the qlist of first-coordinates of the first
graphic element of $P$ (which is the actual Polygon; note that $P$ is
a Graphics object, which consists of a single polygon):

    sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))

We combine together different graphics objects using "+":

    sage: H = G + P + Q
    sage: H
    Graphics object consisting of 3 graphics primitives
    sage: type(H)
    <class 'sage.plot.plot.Graphics'>
    sage: H[1]
    Polygon defined by 3 points
    sage: list(H[1])
    [(1.0, 2.0), (5.0, 6.0), (5.0, 0.0)]

We can put text in a graph:

    sage: L = [[cos(pi*i/100)^3,sin(pi*i/100)] for i in range(200)]
    sage: p = line(L, rgbcolor=(1/4,1/8,3/4))
    sage: t = text('A Bulb', (1.5, 0.25))
    sage: x = text('x axis', (1.5,-0.2))
    sage: y = text('y axis', (0.4,0.9))
    sage: g = p+t+x+y
    sage.: g.show(xmin=-1.5, xmax=2, ymin=-1, ymax=1)

We plot the Riemann zeta function along the critical line and
see the first few zeros:

    sage: p1 = plot(lambda t: arg(zeta(0.5+t*I)), 1,27,rgbcolor=(0.8,0,0))
    sage: p2 = plot(lambda t: abs(zeta(0.5+t*I)), 1,27,rgbcolor=hue(0.7))
    sage: p1+p2
    Graphics object consisting of 2 graphics primitives

Here is a pretty graph:
    sage: g = Graphics()
    sage: for i in range(60):
    ...    p = polygon([(i*cos(i),i*sin(i)), (0,i), (i,0)],\
    ...                rgbcolor=hue(i/40+0.4), alpha=0.2)
    ...    g = g + p
    ...
    sage.: g.show(dpi=200, axes=False)

AUTHORS:
    -- Alex Clemesha and William Stein (2006-04-10): initial version
    -- David Joyner: examples
    -- Alex Clemesha (2006-05-04) major update
    -- Willaim Stein (2006-05-29): fine tuning, bug fixes, better server integration
    -- William Stein (2006-07-01): misc polish
    -- Alex Clemesha (2006-09-29): added contour_plot, frame axes, misc polishing
    -- Robert Miller (2006-10-30): tuning, NetworkX primitive

TODO:
    [] more arithmetic operations on plot objects, e.g., rescaling,
       negation, etc.
    [] ability to change all properties of a graph easily, e.g.,
       the rgbcolor.
"""

#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
#
# WARNING:  Four lines in matplotlib's axes.py module have been
# commented out to make the Point axis work right.
# See the 'scatter' function (page 65) in matplotlib/axes.py
# There are four lines already there that explain the origin
# of the problem, put there buy the author of the function.
# Comment out the *next four lines (approx. 3101-3104)* right below
# the explanation of the hack that produces the problem.
#
# Also:
# On line 3053 of axes.py in the "scatter" function there is a change:
# "if not" -> "if False and not"
#
# NOTE: All this makes me (Alex) wonder if there is a better way to
# make points, but the scatter function provides some
# 'convenient' functionality, so I guess it is worth it?
#
# This is for version 0.86 matplotlib, and these
# changes are included by default in the matplotlib
# included with SAGE.
#
#*****************************************************************************

__doc_exclude = ['SageObject', 'hsv_to_rgb', 'FigureCanvasAgg',\
                 'Figure', 'patches', \
                 'find_axes', 'to_float_list']

DEFAULT_FIGSIZE=[5,4]
DEFAULT_DPI = 125
EMBEDDED_MODE = False
SHOW_DEFAULT = False

do_verify = True

from sage.structure.sage_object import SageObject
import sage.misc.viewer
import sage.misc.misc
verbose = sage.misc.misc.verbose
xsrange = sage.misc.misc.xsrange

import random #for plot adaptive refinement
import os #for viewing and writing images
from colorsys import hsv_to_rgb #for the hue function
from math import sin, cos, modf,pi #for hue and polar_plot

############### WARNING ###
# Try not to import any matplotlib stuff here -- matplotlib is
# slow to import.  (I did benchmarking and found that by not
# importing here, and instead importing when needed below, that
# SAGE startup times are much improved.)  - William
###############

import matplotlib.patches as patches
import matplotlib.cm
from axes import find_axes

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
        sage: G = Graphics(); G
        Graphics object consisting of 0 graphics primitives
        sage: c = circle((1,1), 1)
        sage: G+=c; G
        Graphics object consisting of 1 graphics primitive

    Here we make a graphic of embeded isoceles triangles,
    coloring each one with a different color as we go:

        sage: h=10; c=0.4; p=0.1;
        sage: G = Graphics()
        sage: for x in srange(1,h+1):
        ...        l = [[0,x*sqrt(3)],[-x/2,-x*sqrt(3)/2],[x/2,-x*sqrt(3)/2],[0,x*sqrt(3)]]
        ...        G+=line(l,rgbcolor=hue(c + p*(x/h)))
        sage.: G.show(figsize=[5,5])

    """

    def __init__(self):
        self.__xmin = -1
        self.__xmax = 1
        self.__ymin = -1
        self.__ymax = 1
        self.__fontsize = 6
        self.__show_axes = True
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

    def xmax(self, new=None):
        """
        sage: G = Graphics(); G
        Graphics object consisting of 0 graphics primitives
        sage: G.xmax()
        1
        """
        if new is None:
            return self.__xmax
        self.__xmax = new

    def xmin(self, new=None):
        """
        sage: G = Graphics(); G
        Graphics object consisting of 0 graphics primitives
        sage: G.xmin()
        -1
        """
        if new is None:
            return self.__xmin
        self.__xmin = new

    def ymax(self, new=None):
        """
        sage: G = Graphics(); G
        Graphics object consisting of 0 graphics primitives
        sage: G.ymax()
        1
        """
        if new is None:
            return self.__ymax
        self.__ymax = new

    def ymin(self, new=None):
        """
        sage: G = Graphics(); G
        Graphics object consisting of 0 graphics primitives
        sage: G.ymin()
        -1
        """
        if new is None:
            return self.__ymin
        self.__ymin = new

    def _repr_(self):
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
            sage: G = circle((1,1),2) + circle((2,2),5); G
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

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); G
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
            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); G
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
            sage: del(G[2])
            sage: G
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
            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); G
            Graphics object consisting of 3 graphics primitives

            sage: p = polygon([[1,3],[2,-2],[1,1],[1,3]]);p
            Graphics object consisting of 1 graphics primitive

            sage: G[1] = p[0];G
            Graphics object consisting of 3 graphics primitives

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
            sage: h1 = lambda x : sqrt(x^3  - 1)
            sage: h2 = lambda x : -sqrt(x^3  - 1)
            sage: g1 = plot(h1, 1, 5)
            sage: g2 = plot(h2, 1, 5)
            sage: h = g1 + g2
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

    def append(self, primitive):
        """
        Append an arbitrary GraphicPrimitive to a Graphics object
        """
        if not isinstance(primitive, GraphicPrimitive):
            raise TypeError, "primitive (=%s) must be a GraphicPrimitive"%primitive
        self.__objects.append(primitive)

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

    def _add_xy_axes(self, subplot, xmin, xmax, ymin, ymax, axes_label=None):
        """
        \code{_add_xy_axes} is used when the 'save' method
        of any Graphics object is called.

        additionally this function uses the function 'find_axes'
        from axis.py which attempts to find aesthetically pleasing
        tick and label spacing values.

        some definitons of constants:

        y_axis_xpos : "where on the x-axis to draw the y-axis"
        xstep : "the spacing between major tick marks"
        xtslminor : "x-axis minor tick step list"
        xtslmajor : "x-axis major tick step list"
        yltheight : "where the top of the major ticks go"
        ystheight : "where the top of the minor ticks go"
        ylabel : "where the ylabel is drawn"
        xlabel : "where the xlabel is drawn"

        """
        xmin = float(xmin); xmax=float(xmax); ymin=float(ymin); ymax=float(ymax)
        yspan = ymax - ymin
        xspan = xmax - xmin

        #evalute find_axes for x values and y ticks
        y_axis_xpos, xstep, xtslminor, xtslmajor = find_axes(xmin, xmax)
        yltheight = 0.015 * xspan
        ystheight = 0.25  * yltheight
        ylabel    = y_axis_xpos -2*ystheight

        #evalute find_axes for y values and x ticks
        x_axis_ypos, ystep, ytslminor, ytslmajor = find_axes(ymin, ymax)
        xltheight = 0.015 * yspan
        xstheight = 0.25  * xltheight
        xlabel    = x_axis_ypos - xltheight

        #the x axis line
        subplot.add_line(patches.Line2D([xmin, xmax], [x_axis_ypos, x_axis_ypos],
                                        color='k', linewidth=0.6))

        #the y axis line
        subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos],[ymin, ymax],
                                        color='k', linewidth=0.6))

        def format(z):
            s = str(z)
            if s[-2:] == '.0':
                return s[:-2]
            return s

        #the x-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for x in xtslmajor:
            if x == y_axis_xpos:
                continue
            if self.fontsize() > 0:
                subplot.text(x, xlabel, format(x), fontsize=self.fontsize(),
                             horizontalalignment="center", verticalalignment="top")
            subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xltheight],
                        color='k',linewidth=0.6))

        #now draw the x-axis minor tick marks
        for x in xtslminor:
            subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xstheight],
                        color='k', linewidth=0.6))


        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for y in ytslmajor:
            if y == x_axis_ypos:
                continue
            if self.fontsize() > 0:
                subplot.text(ylabel, y, format(y), fontsize=self.fontsize(), verticalalignment="center",
                        horizontalalignment="right")
            subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + yltheight], [y, y],
                    color='k', linewidth=0.6))

        #now draw the x-axis minor tick marks
        for y in ytslminor:
            subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + ystheight], [y, y],
                color='k',linewidth=0.6))

        #now draw the x and y axis labels
        if axes_label:
            al = axes_label
            if not isinstance(al, (list,tuple)) or len(al) != 2:
                raise TypeError, "axes_label must be a list of two strings."
            #x-axis label
            subplot.text(xmax + 0.2*xstep, x_axis_ypos, str(al[0]), fontsize=6,
                         horizontalalignment="center", verticalalignment="center")
            #y-axis label
            subplot.text(y_axis_xpos, ymax + 0.2*ystep, str(al[1]), fontsize=6,
                         horizontalalignment="center", verticalalignment="center")

    def _add_xy_frame_axes(self, subplot, xmin, xmax, ymin, ymax,
                            axes_with_no_ticks=False, axes_label=None):
        """
        This function is similiar to the above _add_xy_axes but
        it adds a frame with ticks and tick values.
        """
        xmin = float(xmin); xmax=float(xmax); ymin=float(ymin); ymax=float(ymax)
        yspan = ymax - ymin
        xspan = xmax - xmin

        #evalute find_axes for x values and y ticks
        y_axis_xpos, xstep, xtslminor, xtslmajor = find_axes(xmin, xmax)
        #y_axis_xpos = 0 #set to zero for now
        yltheight = 0.015 * xspan
        ystheight = 0.25  * yltheight
        #ylabel    = y_axis_xpos - 2*ystheight
        ylabel    = -2*ystheight

        #evalute find_axes for y values and x ticks
        x_axis_ypos, ystep, ytslminor, ytslmajor = find_axes(ymin, ymax)
        #x_axis_ypos = 0 #set to zero for now
        xltheight = 0.015 * yspan
        xstheight = 0.25  * xltheight
        #xlabel    = x_axis_ypos - xltheight
        xlabel    = -xltheight

        #scale the axes out from the actual plot
        ys = 0.02*yspan
        xs = 0.02*xspan
        ymins = ymin - ys
        ymaxs = ymax + ys
        xmins = xmin - xs
        xmaxs = xmax + xs

        #border horizontal axis:
        #bottom:
        subplot.add_line(patches.Line2D([xmins, xmaxs], [ymins, ymins],
                                        color='k', linewidth=0.6))

        #top:
        subplot.add_line(patches.Line2D([xmins, xmaxs], [ymaxs, ymaxs],
                                        color='k', linewidth=0.6))

        #border vertical axis:
        #left:
        subplot.add_line(patches.Line2D([xmins, xmins], [ymins, ymaxs],
                                        color='k', linewidth=0.6))

        #right:
        subplot.add_line(patches.Line2D([xmaxs, xmaxs], [ymins, ymaxs],
                                        color='k', linewidth=0.6))

        if axes_with_no_ticks:
            #the x axis line
            subplot.add_line(patches.Line2D([xmins, xmaxs], [x_axis_ypos, x_axis_ypos],
                                            color='k', linewidth=0.6))

            #the y axis line
            subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos],[ymins, ymaxs],
                                            color='k', linewidth=0.6))
        def format(z):
            s = str(z)
            if s[-2:] == '.0':
                return s[:-2]
            return s

        #the x-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for x in xtslmajor:
            subplot.text(x, xlabel + ymins, format(x), fontsize=5,
                         horizontalalignment="center", verticalalignment="top")

        #now draw the x-axis minor tick marks
        for x in xtslminor:
            subplot.add_line(patches.Line2D([x, x], [ymins, xstheight + ymins],
                        color='k', linewidth=0.6))
            subplot.add_line(patches.Line2D([x, x], [ymaxs, ymaxs - xstheight],
                        color='k', linewidth=0.6))


        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for y in ytslmajor:
            subplot.text(ylabel + xmins, y, format(y), fontsize=5, verticalalignment="center",
                    horizontalalignment="right")

        #now draw the x-axis minor tick marks
        for y in ytslminor:
            subplot.add_line(patches.Line2D([xmins, ystheight + xmins], [y, y],
                color='k',linewidth=0.6))
            subplot.add_line(patches.Line2D([xmaxs, xmaxs - ystheight], [y, y],
                color='k',linewidth=0.6))


    def _plot_(self, **args):
        return self

    def show(self, xmin=None, xmax=None, ymin=None, ymax=None,
             figsize=DEFAULT_FIGSIZE, filename=None,
             dpi=DEFAULT_DPI, axes=None, axes_label=None,frame=False,
             fontsize=None,
             **args):
        """
        Show this graphics image with the default image viewer.

        EXAMPLES:
            sage: c = circle((1,1), 1, rgbcolor=(1,0,0))
            sage.: c.show(xmin=-1, xmax=3, ymin=-1, ymax=3)

        To correct the apect ratio of certain graphics, it is necessary
        to show with a 'figsize' of square dimensions.

            sage.: c.show(figsize=[5,5], xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can either turn off the drawing of the axes:

            sage.: show(plot(sin,-4,4), axes=False)

        Or you can turn on the drawing of a frame around the plots:

            sage.: show(plot(sin,-4,4), frame=True)

        """
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

    def save(self, filename=None, xmin=None, xmax=None,
             ymin=None, ymax=None, figsize=DEFAULT_FIGSIZE,
             fig=None, sub=None, savenow=True, dpi=DEFAULT_DPI,
             axes=None, axes_label=None, fontsize=None,
             frame=False, verify=True):
        """
        Save the graphics to an image file of type: PNG, PS, EPS, SVG, SOBJ,
        depending on the file extension you give the filename.
            Extension types can be: '.png', '.ps', '.eps', '.svg',
            and '.sobj' (for a SAGE object you can load later).

        EXAMPLES:
            sage: c = circle((1,1),1,rgbcolor=(1,0,0))
            sage.: c.show(xmin=-1,xmax=3,ymin=-1,ymax=3)

            To correct the apect ratio of certain graphics, it is necessary
            to show with a 'figsize' of square dimensions.

            sage.: c.show(figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)
	"""
        global do_verify
        do_verify = verify

        from matplotlib.figure import Figure
        if filename is None:
            i = 0
            while os.path.exists('sage%s.png'%i):
                i += 1
            filename = 'sage%s.png'%i

        try:
            ext = os.path.splitext(filename)[1].lower()
        except IndexError:
            raise ValueError, "file extension must be either 'png', 'eps', 'svg' or 'sobj'"

        if ext == '' or ext == '.sobj':
            SageObject.save(self, filename)
            return

        self.fontsize(fontsize)

        figure = fig
        if not figure:
            figure = Figure(figsize)

        #The line below takes away the excessive whitespace around
        #images.  ('figsize' and  'dpi' still work as expected):
        figure.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)

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
        for g in self.__objects:
            if isinstance(g, GraphicPrimitive_ContourPlot):
                contour = True
            g._render_on_subplot(subplot)

        #adjust the xy limits and draw the axes:
        if axes is None:
            axes = self.__show_axes
        if not contour and axes:
            if frame: #add the frame axes with centered axes with no ticks
                xmin, xmax = self.__xmin, self.__xmax
                ymin, ymax = self.__ymin, self.__ymax
                subplot.set_xlim([xmin - 0.05*abs(xmax - xmin), xmax + 0.05*abs(xmax - xmin)])
                subplot.set_ylim([ymin - 0.05*abs(ymax - ymin), ymax + 0.05*abs(ymax - ymin)])
                self._add_xy_frame_axes(subplot, xmin, xmax, ymin, ymax,
                                        axes_with_no_ticks=True, axes_label=axes_label)
            else:
                xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
                subplot.set_xlim(xmin, xmax)
                subplot.set_ylim(ymin, ymax)
                self._add_xy_axes(subplot, xmin, xmax, ymin, ymax, axes_label=axes_label)

        elif contour:
            xmin, xmax = self.__xmin, self.__xmax
            ymin, ymax = self.__ymin, self.__ymax
            subplot.set_xlim([xmin - 0.05*abs(xmax - xmin), xmax + 0.05*abs(xmax - xmin)])
            subplot.set_ylim([ymin - 0.05*abs(ymax - ymin), ymax + 0.05*abs(ymax - ymin)])
            if axes:
                self._add_xy_frame_axes(subplot, xmin, xmax, ymin, ymax, axes_label=axes_label)
        else:
            xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)
            subplot.set_xlim(xmin, xmax)
            subplot.set_ylim(ymin, ymax)


        # you can output in PNG, PS, or SVG format, depending on the file extension
        if savenow:
            if ext in ['.eps', '.ps']:
                from matplotlib.backends.backend_ps import FigureCanvasPS
                canvas = FigureCanvasPS(figure)
                if dpi is None:
                    dpi = 72
            elif ext == '.svg':
                from matplotlib.backends.backend_svg import FigureCanvasSVG
                if dpi is None:
                    dpi = 80
                canvas = FigureCanvasSVG(figure)
            elif ext == '.png':
                from matplotlib.backends.backend_agg import FigureCanvasAgg
                canvas = FigureCanvasAgg(figure)
                if dpi is None:
                    dpi = 150
            else:
                raise ValueError, "file extension must be either 'png', 'eps', 'svg' or 'sobj'"
            #canvas.resize(100,100)
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
                    verbose("WARNING: Ignoring option '%s'=%s"%(k,O[k]), level=0)
                    t = True
            if t:
                s = "\nThe allowed options for %s are:\n"%self
                K.sort()
                for k in K:
                    if A.has_key(k):
                        s += "    %-15s%-60s\n"%(k,A[k])
                verbose(s, level=0)


        if 'hue' in O:
            t = O['hue']
            if not isinstance(t, (tuple,list)):
                t = [t,1,1]
            O['rgbcolor'] = hue(*t)
            del O['hue']
        return O

    def _repr_(self):
        return "Graphics primitive"



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
                'hue':'The color given as a hue.'}

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
        options = self.options()
        p = patches.Line2D(self.xdata, self.ydata)
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
                'resolution': '',
                'fill': 'Whether or not to fill the polygon.',
                'thickness':'How thick the border of the polygon is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.'}

    def _repr_(self):
        return "Circle defined by (%s,%s) with r=%s"%(self.x, self.y, self.r)

    def _render_on_subplot(self, subplot):
        options = self.options()
        res = int(options['resolution'])
        p = patches.Circle((float(self.x), float(self.y)), float(self.r),resolution=res)
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
                'fill':'Fill contours or not'}

    def _repr_(self):
        return "ContourPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        options = self.options()
        fill = options['fill']
        cmap = options['cmap']
        try:
            cmap = matplotlib.cm.__dict__[cmap]
        except KeyError:
            from matplotlib.colors import LinearSegmentedColormap as C
            possibilities = ', '.join([str(x) for x in matplotlib.cm.__dict__.keys() if \
                                       isinstance(matplotlib.cm.__dict__[x], C)])
            print "The possible color maps include: %s"%possibilities
            raise RuntimeError, "Color map %s not known"%cmap

        x0,x1 = float(self.xrange[0]), float(self.xrange[1])
        y0,y1 = float(self.yrange[0]), float(self.yrange[1])
        if fill:
            subplot.contourf(self.xy_data_array, cmap=cmap, extent=(x0,x1,y0,y1))
        else:
            subplot.contour(self.xy_data_array, cmap=cmap, extent=(x0,x1,y0,y1))

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
        options = self.options()
        deg1 = self.rad1*(360.0/(2.0*pi))
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
        #see top of this file for Point info
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
        if i == 0:
            return self.xdata
        elif i == 1:
            return self.ydata
        else:
            raise IndexError, "Index out of range"

    def _render_on_subplot(self,subplot):
        options = self.options()
        c = to_mpl_color(options['rgbcolor'])
        a = float(options['alpha'])
        s = int(options['pointsize'])
        faceted = options['faceted'] #faceted=True colors the edge of point
        subplot.scatter(self.xdata, self.ydata, s, c, alpha=a, faceted=False)

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
                'vertical_alignment':'if True align vertically.',
                'horizontal_alignment':'if True align vertically.'}

    def _render_on_subplot(self, subplot):
        options = self.options()
        c = to_mpl_color(options['rgbcolor'])
        f = int(options['fontsize'])
        va = options['vertical_alignment']
        ha = options['horizontal_alignment']
        subplot.text(float(self.x), float(self.y), self.string, color=c, fontsize=f,
                        verticalalignment=va,horizontalalignment=ha)

class GraphicPrimitive_NetworkXGraph(GraphicPrimitive):
    """
    Primitive class that initializes the NetworkX graph type

    INPUT:
        graph -- a NetworkX graph
        pos -- an optional positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {   0: [-0.91679746, 0.88169588,],
                1: [ 0.47294849, 1.125     ,],
                2: [ 1.125     ,-0.12867615,],
                3: [ 0.12743933,-1.125     ,],
                4: [-1.125     ,-0.50118505,]   }
        with_labels -- determines whether labels for nodes are plotted
        node_size -- node size

    EXAMPLE:
        sage: import networkx as NX
        sage: D = NX.dodecahedral_graph()
        sage: NGP = GraphicPrimitive_NetworkXGraph(D)
        sage: g = Graphics()
        sage: g.append(NGP)
        sage: g.axes(False)
        sage: g.save('a.png')
    """
    def __init__(self, graph, pos=None, with_labels=True, node_size=300):
        self.__nxg = graph
        if len(self.__nxg) != 0:
            import networkx as NX
            self.__node_size = node_size
            self.__with_labels = with_labels
            if pos is None:
                self.__pos = NX.drawing.spring_layout(self.__nxg)
            else:
                self.__pos = pos

            # adjust the plot
            # TODO: right now, the bounds are assumed to be +-1
            nodelist=self.__nxg.nodes()
            xes = [self.__pos[v][int(0)] for v in nodelist]
            ys = [self.__pos[v][int(1)] for v in nodelist]
            xmin = min(xes)
            xmax = max(xes)
            ymin = min(ys)
            ymax = max(ys)
            st = float(0.125) # scaling term: looks better, maybe should
                              # depend on node_size?
            if xmax == xmin:
                xmax += 1
                xmin -= 1
            if ymax == ymin:
                ymax += 1
                ymin -= 1
            for v in nodelist:
                self.__pos[v][int(0)] = ((2 + (2*st))/(xmax-xmin))*(self.__pos[v][int(0)] - xmax) + st + 1
                self.__pos[v][int(1)] = ((2 + (2*st))/(ymax-ymin))*(self.__pos[v][int(1)] - ymax) + st + 1

    def _render_on_subplot(self, subplot):
        if len(self.__nxg) != 0:
            import networkx as NX
            node_size = float(self.__node_size)
            NX.draw_networkx_nodes(G=self.__nxg, pos=self.__pos, ax=subplot, node_size=node_size)
            NX.draw_networkx_edges(G=self.__nxg, pos=self.__pos, ax=subplot, node_size=node_size)
            if self.__with_labels:
                labels = {}
                for v in self.__nxg:
                    labels[v] = str(v)
                NX.draw_networkx_labels(self.__nxg, self.__pos, labels=labels, ax=subplot)

######################################################################
#                                                                    #
#    Graphics Primitives Factories -- construct GraphicPrimitives    #
#                                                                    #
######################################################################

class GraphicPrimitiveFactory:
    def __init__(self):
        # options for this specific graphics primitive.
        self.reset()

    def reset(self):
        # First the default options for all graphics primitives
        self.options = {'alpha':1,'thickness':1,'rgbcolor':(0,0,0)}
        self._reset()

    def _coerce(self, xdata, ydata):
        return to_float_list(xdata), to_float_list(ydata)

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
        #when contour_plot is called??  Here we go:
        plot_points = int(options['plot_points'])
        xstep = abs(float(xrange[0]) - float(xrange[1]))/plot_points
        ystep = abs(float(yrange[0]) - float(yrange[1]))/plot_points
        xy_data_array = [[float(f(x, y)) for x in xsrange(xrange[0], xrange[1], xstep)]
                                     for y in xsrange(yrange[0], yrange[1], ystep)]
        return self._from_xdata_ydata(xy_data_array, xrange, yrange, options=options)

class GraphicPrimitiveFactory_disk(GraphicPrimitiveFactory):
    def __call__(self, point, radius, angle, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata((float(point[0]), float(point[1])),float(radius),
                    (float(angle[0]), float(angle[1])), options=options)

class GraphicPrimitiveFactory_text(GraphicPrimitiveFactory):
    def __call__(self, string, point, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata(string, (float(point[0]), float(point[1])), options=options)

class GraphicPrimitiveFactory_from_point_list(GraphicPrimitiveFactory):
    def __call__(self, points, coerce=True, **kwds):
        try:
            return points._plot_(**kwds)
        except AttributeError:
            pass
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v

        done = False
        if not isinstance(points, (list,tuple)) or \
           (isinstance(points,(list,tuple)) and len(points) == 2):
            try:
                xdata = [float(points[0])]
                ydata = [float(points[1])]
                done = True
            except TypeError:
                pass

        if not done:
            if coerce:
                xdata = []
                ydata = []
                for z in points:
                    xdata.append(float(z[0]))
                    ydata.append(float(z[1]))
            else:
                xdata = [z[0] for z in points]
                ydata = [z[1] for z in points]

        return self._from_xdata_ydata(xdata, ydata, True, options=options)

class LineFactory(GraphicPrimitiveFactory_from_point_list):
    r"""
    Create the line through the given list of points.

    Type line.options for a dictionary of the default options for
    lines.  You can change this to change the defaults for all future
    lines.  Use line.reset() to reset to the default options.

    EXAMPLES:
    A blue conchoid of Nicomedes:

        sage: L = [[1+5*cos(pi/2+pi*i/100), tan(pi/2+pi*i/100)*(1+5*cos(pi/2+pi*i/100))] for i in range(1,100)]
        sage: p = line(L, rgbcolor=(1/4,1/8,3/4))

    A blue hypotrochoid (3 leaves):

        sage: n = 4; h = 3; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: p = line(L, rgbcolor=(1/4,1/4,3/4))

    A blue hypotrochoid (4 leaves):

        sage: n = 6; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: p = line(L, rgbcolor=(1/4,1/4,3/4))

    A red limacon of Pascal:

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,101)]
        sage: p = line(L, rgbcolor=(1,1/4,1/2))

    A light green trisectrix of Maclaurin:

        sage: L = [[2*(1-4*cos(-pi/2+pi*i/100)^2),10*tan(-pi/2+pi*i/100)*(1-4*cos(-pi/2+pi*i/100)^2)] for i in range(1,100)]
        sage: p = line(L, rgbcolor=(1/4,1,1/8))

    A green lemniscate of Bernoulli:

        sage: v = [(1/cos(-pi/2+pi*i/100), tan(-pi/2+pi*i/100)) for i in range(201)]
        sage: L = [(a/(a^2+b^2), b/(a^2+b^2)) for a,b in v]
        sage: p = line(L, rgbcolor=(1/4,3/4,1/8))

    A red plot of the Jacobi elliptic function $\text{sn}(x,2)$, $-3<x<3$:

        sage: L = [(i/100.0, maxima.eval('jacobi_sn (%s/100.0,2.0)'%i)) for i in range(-300,300)]
        sage: p = line(L, rgbcolor=(3/4,1/4,1/8))

    A red plot of $J$-Bessel function $J_2(x)$, $0<x<10$:

        sage: L = [(i/10.0, maxima.eval('bessel_j (2,%s/10.0)'%i)) for i in range(100)]
        sage: p = line(L, rgbcolor=(3/4,1/4,5/8))

    A purple plot of the Riemann zeta function $\zeta(1/2 + it)$, $0<t<30$:

        sage: v = [zeta(0.5 + i/10 * I) for i in range(300)]
        sage: L = [(z.real(), z.imag()) for z in v]
        sage: p = line(L, rgbcolor=(3/4,1/2,5/8))

    A purple plot of the Hasse-Weil $L$-function $L(E, 1 + it)$, $-1<t<10$:

        sage: E = EllipticCurve('37a')
        sage: vals = E.Lseries_values_along_line(1-I, 1+10*I, 100) # critical line
        sage: L = [(z[1].real(), z[1].imag()) for z in vals]
        sage: p = line(L, rgbcolor=(3/4,1/2,5/8))

    A red, blue, and green "cool cat":

        sage: ncos = lambda x: -cos(x)
        sage: G = plot(ncos, -2, 2, thickness=5, rgbcolor=(0.5,1,0.5))
        sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))
        sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
        sage: H = G + P + Q
    """
    def _reset(self):
        self.options = {'alpha':1,'rgbcolor':(0,0,0),'thickness':1}

    def _repr_(self):
        return "type line? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g.append(GraphicPrimitive_Line(xdata, ydata, options))
        try:
            g._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass
        return g

# unique line instance
line = LineFactory()


class CircleFactory(GraphicPrimitiveFactory_circle):
    """
    A circle at a point = (x,y) with radius = r
    Type circle.options to see all options

    EXAMPLES:
    sage: c = circle((1,1),1,rgbcolor=(1,0,0))
    sage.: c.show(xmin=-1,xmax=3,ymin=-1,ymax=3)

    To correct the apect ratio of certain graphics, it is necessary
    to show with a 'figsize' of square dimensions.

    sage.: c.show(figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)

    Here we make an more complicated plot with many circles of different colors

    sage: g = Graphics()
    sage: step=6; ocur=1/5; paths=16;
    sage: for r in range(1,paths+1):
    ...       for x,y in [((r+ocur)*cos(n), (r+ocur)*sin(n)) for n in srange(0, 2*pi+pi/step, pi/step)]:
    ...           g += circle((x,y), ocur, rgbcolor=hue(r/paths))
    ...       rnext = (r+1)^2
    ...       ocur = (rnext-r)-ocur
    ...
    sage.: g.show(xmin=-(paths+1)^2, xmax=(paths+1)^2, ymin=-(paths+1)^2, ymax=(paths+1)^2, figsize=[6,6])
    """
    def _reset(self):
        self.options={'alpha':1,'fill':False,'thickness':1,'rgbcolor':(0, 0, 0),'resolution':40}

    def _repr_(self):
        return "type circle? for help and examples"

    def _from_xdata_ydata(self, point, r, options):
        g = Graphics()
        x = float(point[0])
        y = float(point[1])
        r = float(r)
        g.append(GraphicPrimitive_Circle(x, y, r, options))
        g._extend_axes(x+r, x-r, y+r, y-r)
        return g


#an unique circle instance
circle = CircleFactory()

class ContourPlotFactory(GraphicPrimitiveFactory_contour_plot):
    """
    \code{contour_plot} takes a function of two variables, f(x,y)
    and plots contour lines of the function over the specified
    xrange and yrange as demonstrated below.
    contour_plot(f, (xmin, xmax), (ymin, ymax))

    EXAMPLES:

        Here we plot a simple funtion of two variables:
        sage: def f(x,y):
        ...       return cos(x^2 + y^2)
        sage: contour_plot(f, (-4, 4), (-4, 4))
         Graphics object consisting of 1 graphics primitive


        Here we change the ranges and add some options:
        sage: def h(x,y):
        ...       return (x^2)*cos(x*y)
        sage: contour_plot(h, (-10, 5), (-5, 5), fill=False, plot_points=100)
         Graphics object consisting of 1 graphics primitive

    """
    def _reset(self):
        self.options={'plot_points':50, 'fill':True, 'cmap':'gray'}

    def _repr_(self):
        return "type contour_plot? for help and examples"

    def _from_xdata_ydata(self, xy_data_array, xrange, yrange, options):
        g = Graphics()
        g.__xmin = xrange[0]
        g.__xmax = xrange[1]
        g.__ymin = yrange[0]
        g.__ymax = yrange[1]
        g.append(GraphicPrimitive_ContourPlot(xy_data_array, xrange, yrange, options))
        g._extend_axes(xrange[0], xrange[1], yrange[0], yrange[1])
        return g

#unique contour_plot instance
contour_plot = ContourPlotFactory()

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
        sage.: P.show(figsize=(4,4),xmin=-2,xmax=2,ymin=-2,ymax=2)

    """
    def _reset(self):
        self.options={'alpha':1,'fill':True,'rgbcolor':(0,0,0),'thickness':0}

    def _repr_(self):
        return "type disk? for help and examples"

    def _from_xdata_ydata(self, point, r, angle, options):
        g = Graphics()
        g.append(GraphicPrimitive_Disk(point, r, angle, options))
        g._extend_axes(2*r, -2*r, 2*r, -2*r)
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
        sage: p1 = point((0.5,0.5), rgbcolor=hue(0.75))

        Here are some random larger red points, given as a list of tuples
        sage: p2 = point(((0.5,0.5),(1,2),(0.5,0.9),(-1,-1)),rgbcolor=hue(1),pointsize=30)

    """
    def _reset(self):
        self.options = {'alpha':1,'pointsize':10,'faceted':False,'rgbcolor':(0,0,0)}

    def _repr_(self):
        return "type point? for options help"

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g.append(GraphicPrimitive_Point(xdata, ydata, options))
        try:
            g._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))
        except ValueError:
            pass
        return g

# unique point instance
point = PointFactory()

class TextFactory(GraphicPrimitiveFactory_text):
    """
    Text at the point (x,y)
    Type text.options for a dictionary of options.

    EXAMPLES:
        Type this to see some text in top right plane:

        sage: t = text("SAGE is really neat!!",(2,12))

        Type this to see the same text in larger font and colored red:

        sage: t = text("SAGE is really neat!!",(2,12),fontsize=20,rgbcolor=(1,0,0))

        You can also center text differently:

        sage: t1 = text("Hello",(1,1), vertical_alignment="top")
        sage: t2 = text("World", (1,0.5), horizontal_alignment="left")

    """
    def _reset(self):
        self.options = {'fontsize':10, 'rgbcolor':(0,0,0),
                        'horizontal_alignment':'center',
                        'vertical_alignment':'center'}

    def _repr_(self):
        return "type text? for help and examples"

    def _from_xdata_ydata(self, string, point, options):
        g = Graphics()
        g.append(GraphicPrimitive_Text(string, point, options))
        g._extend_x_axis(2*point[0])
        g._extend_y_axis(2*point[1])
        return g

# unique text instance
text = TextFactory()

class PolygonFactory(GraphicPrimitiveFactory_from_point_list):
    """
    Type polygon.options for a dictionary of the default
    options for polygons.  You can change this to change
    the defaults for all future polygons.  Use polygon.reset()
    to reset to the default options.

    EXAMPLES:
    We create a purple-ish polygon:
        sage: polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))
        Graphics object consisting of 1 graphics primitive

    Some modern art -- a random polygon:
        sage: v = [(randrange(-5,5), randrange(-5,5)) for _ in range(10)]
        sage: polygon(v)
        Graphics object consisting of 1 graphics primitive

    A purple hexagon:

        sage: L = [[cos(pi*i/3),sin(pi*i/3)] for i in range(6)]
        sage: p = polygon(L, rgbcolor=(1,0,1))

    A green deltoid:

        sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))

    A blue hypotrochoid:

        sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))

    Another one:

        sage: n = 4; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/8,1/4,3/4))

    A purple epicycloid:

        sage: m = 9; b = 1
        sage: L = [[m*cos(pi*i/100)+b*cos((m/b)*pi*i/100),m*sin(pi*i/100)-b*sin((m/b)*pi*i/100)] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(7/8,1/4,3/4))

    A brown astroid:

        sage: L = [[cos(pi*i/100)^3,sin(pi*i/100)^3] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(3/4,1/4,1/4))

    And, my favorite, a greenish blob:

        sage: L = [[cos(pi*i/100)*(1+cos(pi*i/50)), sin(pi*i/100)*(1+sin(pi*i/50))] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/8, 3/4, 1/2))

    This one is for my wife:

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,100)]
        sage: p = polygon(L, rgbcolor=(1,1/4,1/2))

    AUTHORS:
        -- David Joyner (2006-04-14): the long list of examples above.
    """
    def _reset(self):
        self.options={'alpha':1,'rgbcolor':(0,0,0),'thickness':0}

    def _repr_(self):
        return "SAGE polygon; type polygon? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g.append(GraphicPrimitive_Polygon(xdata, ydata, options))
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

    where X is a SAGE object that either is callable and returns
    numbers that can be coerced to floats, or has a _plot_ method
    that returns a GraphicPrimitive object.

    Type plot.options for a dictionary of the default
    options for plots.  You can change this to change
    the defaults for all future plots.  Use plot.reset()
    to reset to the default options.

    The options are:
        plot_points -- the number of points to initially plot before
                       doing adaptive refinement
        plot_division -- the maximum number points including those
                       computed during adaptive refinement
        max_bend      -- parameter that affects adaptive refinement

        xmin -- starting x value
        xmax -- ending x value

    Note that this function does NOT simply sample equally spaced
    points between xmin and xmax.  Instead it computes equally spaced
    points and add small perturbations to them.  This reduces the
    possibility of, e.g., sampling sin only at multiples of $2\pi$,
    which would yield a very misleading graph.

    EXAMPLES:
    We plot the sin function:
        sage: P = plot(sin, 0,10); P
        Graphics object consisting of 1 graphics primitive
        sage: len(P)     # number of graphics primitives
        1
        sage: len(P[0])  # how many points were computed
        201
        sage: P = plot(sin, 0,10, plot_points=10); P
        Graphics object consisting of 1 graphics primitive
        sage: len(P[0])   # random output
        80

    Use \code{show(plot(sin, 0,10))} or \code{plot(sin,0,10).show()}
    to show the corresponding graphics object.

    We draw the graph of an elliptic curve as the union
    of graphs of 2 functions.
        sage: def h1(x): return sqrt(x^3  - 1)
        sage: def h2(x): return -sqrt(x^3  - 1)
        sage: plot([h1, h2], 1,4)    # random output because of random sampling
        Graphics object consisting of 2 graphics primitives

    We can also directly plot the elliptic curve:
        sage: E = EllipticCurve([0,-1])
        sage: P = plot(E, 1, 4, rgbcolor=hue(0.6))
    """
    def _reset(self):
        o = self.options
        o['plot_points'] = 200
        o['plot_division'] = 1000      # change when fix adapter
        o['max_bend'] = 0.1            # change when fix adapter

    def _repr_(self):
        return "plot; type plot? for help and examples."

    def __call__(self, funcs, xmin=None, xmax=None, parametric=False,
                 polar=False, show=None, **kwds):
        if show is None:
            show = SHOW_DEFAULT
        try:
            G = funcs._plot_(xmin=xmin, xmax=xmax, **kwds)
            if show:
                G.show(**kwds)
            return G
        except AttributeError:
            pass
        try:
            G = funcs.plot(xmin=xmin, xmax=xmax, **kwds)
            if show:
                G.show(**kwds)
            return G
        except AttributeError:
            pass

        if xmin is None:
            xmin = -1
        if xmax is None:
            xmax = 1  # defaults
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        #check to see if funcs is a list of functions that will
        #be all plotted together.
        if isinstance(funcs, (list, tuple)) and not parametric:
            G = Graphics()
            for i in range(0, len(funcs)):
                G += plot(funcs[i], xmin, xmax, polar=polar, **kwds)
            if show:
                G.show(**kwds)
            return G
        #parametric_plot will be a list or tuple of two functions (f,g)
        #and will plotted as (f(x), g(x)) for all x in the given range
        if parametric:
            f,g = funcs
        #or we have only a single function to be plotted:
        else:
            f = funcs
        xmin = float(xmin)
        xmax = float(xmax)
        plot_points = int(options['plot_points'])
        del options['plot_points']
        delta = (xmax - xmin) / plot_points
        data = []
        dd = delta
        for i in xrange(plot_points + 1):
            x = xmin + i*delta
            if i < plot_points:
                x += delta*random.random()
                if x > xmax:
                    x = xmax
            else:
                x = xmax  # guarantee that we get the last point.
            try:
                y = f(x)
                data.append((x, float(y)))
            except (TypeError, ValueError), msg:
                raise ValueError, "Unable to compute f(%s)"%x
        # adaptive refinement
        i, j = 0, 0
        max_bend = float(options['max_bend'])
        del options['max_bend']
        plot_division = int(options['plot_division'])
        del options['plot_division']
        while i < len(data) - 1:
            if abs(data[i+1][1] - data[i][1]) > max_bend:
                x = (data[i+1][0] + data[i][0])/2
                y = float(f(x))
                data.insert(i+1, (x, y))
                j += 1
                if j > plot_division:
                    # wrong -- just for testing (so no infinite loop)
                    break
            else:
                i += 1
        if parametric:
            data = [(fdata, g(x)) for x, fdata in data]
        if polar:
            data = [(y*cos(x), y*sin(x)) for x, y in data]
        G = line(data, coerce=False, **options)
        if show:
            G.show(**kwds)
        return G

# unique plot instance
plot = PlotFactory()

########## misc functions ###################

#parametric plot
def parametric_plot((f,g), tmin, tmax, show=None, **kwargs):
    """
    parametric_plot takes two functions as a list or a tuple and make
    a plot with the first function giving the x coordinates and the
    second function giving the y coordinates.

    INPUT:
        (f,g) -- tuple of functions
        tmin -- start value of t
        tmax -- end value of t
        show -- whether or not to show the plot immediately (default: True)
        other options -- passed to plot.

    EXAMPLE:
        sage: f = lambda t: sin(t)
        sage: g = lambda t: sin(2*t)
        sage: parametric_plot((f,g),0,2*pi,rgbcolor=hue(0.6))
        Graphics object consisting of 1 graphics primitive
    """
    if show is None:
        show = SHOW_DEFAULT
    return plot((f,g), tmin, tmax, parametric=True, show=show, **kwargs)

#polar plot
def polar_plot(funcs, xmin, xmax, show=None, **kwargs):
    """
    polar_plot takes a single function or a list or tuple of functions
    and plots them parametrically in the given range.

    EXAMPLES:
    Here is a blue 8-leaved petal:
        sage: p1 = polar_plot(lambda x:sin(5*x)^2, 0, 2*pi, rgbcolor=hue(0.6))

    A red figure-8:
        sage: p2 = polar_plot(lambda x:sqrt(1 - sin(x)^2), 0, 2*pi, rgbcolor=hue(1.0))

    A green limacon of Pascal:
        sage: p3 = polar_plot(lambda x:2 + 2*cos(x), 0, 2*pi, rgbcolor=hue(0.3))

    """
    if show is None:
        show = SHOW_DEFAULT
    return plot(funcs, xmin, xmax, polar=True, **kwargs)

#list plot
def list_plot(data, plotjoined=False, show=None, **kwargs):
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
        Graphics object consisting of 1 graphics primitive

        sage: r = [(random(),random()) for _ in range(20)]

    Here are a bunch of random red points:
        sage: list_plot(r,rgbcolor=(1,0,0))
        Graphics object consisting of 1 graphics primitive

    This gives all the random points joined in a purple line:
        sage: list_plot(r, plotjoined=True, rgbcolor=(1,0,1))
        Graphics object consisting of 1 graphics primitive
    """
    if show is None:
        show = SHOW_DEFAULT
    if not isinstance(data[0], (list, tuple)):
        data = zip(range(len(data)),data)
    if plotjoined:
        P = line(data, **kwargs)
    else:
        P = point(data, **kwargs)
    if show:
        P.show(**kwargs)
    return P

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

        sage: p = plot(sin, -2, 2, rgbcolor=hue(0.6))

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

    def _repr_(self):
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

    def __len__(self):
        return len(self._glist)

    def append(self, g):
        self._glist.append(g)

    def _render(self, filename, dpi=None, figsize=DEFAULT_FIGSIZE, **args):
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
            g.save(filename, dpi=dpi, fig=figure,
                   sub=subplot, savenow = (i==dims), verify=do_verify,
                   **args)   # only save if i==dims.

    def save(self, filename=None, dpi=DEFAULT_DPI, figsize=DEFAULT_FIGSIZE, **args):
        """
        save the \code{graphics_array} to
            (for now) a png called 'filename'.
        """
        self._render(filename, dpi=dpi, figsize=figsize, **args)

    def show(self, filename=None, dpi=DEFAULT_DPI, figsize=DEFAULT_FIGSIZE, **args):
        r"""
        Show this graphics array using the default viewer.
        """
        if EMBEDDED_MODE:
            self.save(filename, dpi=dpi, figsize=figsize, **args)
            return
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self._render(filename, dpi=dpi, figsize=figsize, **args)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(
                         sage.misc.viewer.browser(), filename))
        #os.system('gqview %s >/dev/null&'%filename)


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

        sage: f = lambda x: sin(x)
        sage: g = lambda x: sin(2*x)
        sage: h = lambda x: sin(4*x)
        sage: p1 = plot(f,-2*pi,2*pi,rgbcolor=hue(0.5))
        sage: p2 = plot(g,-2*pi,2*pi,rgbcolor=hue(0.9))
        sage: p3 = parametric_plot((f,g),0,2*pi,rgbcolor=hue(0.6))
        sage: p4 = parametric_plot((f,h),0,2*pi,rgbcolor=hue(1.0))

    Now make a graphics array out of the plots;
    Ten you can type either: \code{ga.show()} or \code{ga.save()}.

        sage: ga = graphics_array(((p1,p2),(p3,p4)))

    Here we give only one row:
        sage: p1 = plot(sin,-4,4)
        sage: p2 = plot(cos,-4,4)
        sage: graphics_array([p1, p2])
        Graphics Array of size 1 x 2

    """
    if not n is None:
        # Flatten then reshape input
        n = int(n)
        m = int(m)
        array = reshape(array, n, m)

    G = GraphicsArray(array)
    return G



def arrowhead(x,y,angle=0,spread=0.1,length=0.05,**options):
    """
    Draw an arrowhead with tip at (x,y), rotated the given angle,
    with given spread (in radians) and sides having the given length.

    EXAMPLES:

    """
    s = spread/2; r = length
    v = [(x - r*cos(angle-s), y-r*sin(angle-s)), (x,y), \
         (x - r*cos(angle+s), y-r*sin(angle+s))]
    return line(v, **options)
