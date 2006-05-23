r"""
2D Plotting

SAGE provides 2-d plotting functionality with an interface inspired by
the interface for plotting in Mathematica.  The underlying rendering
is mostly implemented using the matplotlib Python library.

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
    Graphics object consisting of 3 graphics primitives:
            0 -- Line defined by 201 points
            1 -- Polygon defined by 3 points
            2 -- Polygon defined by 3 points

    sage: type(H)
    <class 'sage.plot.plot.Graphics'>
    sage: H[1]
    Polygon defined by 3 points
    sage: list(H[1])
    [(1.0, 2.0), (5.0, 6.0), (5.0, 0.0)]

AUTHORS:
    -- Alex Clemesha and William Stein (2006-04-10): initial version
    -- David Joyner: examples
    -- Alex Clemesha (2006-05-04) major update
"""

__doc_exclude = ['SageObject', 'hsv_to_rgb', 'FigureCanvasAgg', 'Value', \
                 'Figure', 'patches', 'flatten', 'to_float_list']  #no ref manual

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
                 'Value', 'Figure', 'patches', 'flatten', \
                 'find_axes']

DEFAULT_FIGSIZE=[5,4]

import sage.misc.log
import sage.misc.misc

import os

from sage.ext.sage_object import SageObject
from colorsys import hsv_to_rgb #for the hue function
from math import modf

try:

    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.backends.backend_ps import FigureCanvasPS
    from matplotlib.backends.backend_svg import FigureCanvasSVG
    from matplotlib.transforms import Value
    from matplotlib.figure import Figure
    import matplotlib.patches as patches
    from matplotlib.cbook import flatten

except ImportError, msg:
    print msg
    print "WARNING -- matplotlib did not build correctly as part of SAGE."


from axes import find_axes

def is_Graphics(x):
    """
    Return True if x is a Graphics object.

    EXAMPLES:
        sage: is_Graphics(1)
        False
        sage: is_Graphics(disk((0.0,0.0),1,0,90))
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
        Graphics object consisting of 0 graphics primitives:
        sage: c = circle((1,1), 1)
        sage: G+=c; G
        Graphics object consisting of 1 graphics primitives:
                0 -- Circle defined by (1.0,1.0) with r=1.0

        Here we make a graphic of embeded isoceles triangles,
	coloring each one with a different color as we go:

        sage: h=10; c=0.4; p=0.1;
        sage: G = Graphics()
        sage: for x in srange(1,h+1):
        ...        l = [[0,x*sqrt(3)],[-x/2,-x*sqrt(3)/2],[x/2,-x*sqrt(3)/2],[0,x*sqrt(3)]]
        ...        G+=line(l,rgbcolor=hue(c + p*(x/h)))
        sage: G.save("triangle.png",figsize=[5,5])

    """

    def __init__(self):
        self.__xmin = -1
        self.__xmax = 1
        self.__ymin = -1
        self.__ymax = 1
        self.__objects = []

    def xmax(self):
        """
	sage: G = Graphics(); G
	Graphics object consisting of 0 graphics primitives:
	sage: G.xmax()
	1
        """
        return self.__xmax

    def xmin(self):
	"""
	sage: G = Graphics(); G
	Graphics object consisting of 0 graphics primitives:
	sage: G.xmin()
	-1
	"""
        return self.__xmin

    def ymax(self):
	"""
	sage: G = Graphics(); G
	Graphics object consisting of 0 graphics primitives:
	sage: G.ymax()
	1
	"""
        return self.__ymax

    def ymin(self):
	"""
	sage: G = Graphics(); G
	Graphics object consisting of 0 graphics primitives:
	sage: G.ymin()
	-1
	"""
        return self.__ymin

    def _repr_(self):
        pr, i = '', 0
        for x in self:
            pr += '\n\t%s -- %s'%(i, x)
            i += 1
        return "Graphics object consisting of %s graphics primitives:%s"%(
            len(self), pr)

    def __getitem__(self, i):
	"""
	Returns the ith graphics primitive object:

	EXAMPLE:
	    sage: G = circle((1,1),2) + circle((2,2),5); G
	    Graphics object consisting of 2 graphics primitives:
        	    0 -- Circle defined by (1.0,1.0) with r=2.0
        	    1 -- Circle defined by (2.0,2.0) with r=5.0
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
            Graphics object consisting of 3 graphics primitives:
                    0 -- Circle defined by (1.0,1.0) with r=1.0
                    1 -- Circle defined by (1.0,2.0) with r=1.0
                    2 -- Circle defined by (1.0,2.0) with r=5.0
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
            Graphics object consisting of 3 graphics primitives:
                    0 -- Circle defined by (1.0,1.0) with r=1.0
                    1 -- Circle defined by (1.0,2.0) with r=1.0
                    2 -- Circle defined by (1.0,2.0) with r=5.0
            sage: len(G)
            3
            sage: del(G[2])
            sage: G
            Graphics object consisting of 2 graphics primitives:
                    0 -- Circle defined by (1.0,1.0) with r=1.0
                    1 -- Circle defined by (1.0,2.0) with r=1.0
	    sage: len(G)
	    2
        """
        del self.__objects[int(i)]

    def __setitem__(self, i, x):
        """
	You can replace an GraphicsPrimitive (point, line, circle, etc...)
	in a Graphics object G with any other GraphicsPrimitive

	EXAMPLES:
            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); G
            Graphics object consisting of 3 graphics primitives:
                    0 -- Circle defined by (1.0,1.0) with r=1.0
                    1 -- Circle defined by (1.0,2.0) with r=1.0
                    2 -- Circle defined by (1.0,2.0) with r=5.0

	    sage: p = polygon([[1,3],[2,-2],[1,1],[1,3]]);p
	    Graphics object consisting of 1 graphics primitives:
                    0 -- Polygon defined by 4 points

            sage: G[1] = p[0];G
	    Graphics object consisting of 3 graphics primitives:
            	    0 -- Circle defined by (1.0,1.0) with r=1.0
        	    1 -- Polygon defined by 4 points
        	    2 -- Circle defined by (1.0,2.0) with r=5.0

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
            sage: h = g1 + g2; h
	    Graphics object consisting of 2 graphics primitives:
            	    0 -- Line defined by 205 points
        	    1 -- Line defined by 205 points

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

    def _circle(self, x, y, r, options):
        self.__objects.append(GraphicPrimitive_Circle(x, y, r, options))
        self._extend_axes(x+r, x-r, y+r, y-r)

    def _disk(self, point, r, theta1, theta2, options):
        self.__objects.append(GraphicPrimitive_Disk(point, r, theta1, theta2, options))
        self._extend_axes(2*r, -2*r, 2*r, -2*r)

    def _line(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Line(xdata, ydata, options))
        self._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))

    def _point(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Point(xdata, ydata, options))
        self._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))

    def _polygon(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Polygon(xdata, ydata, options))
        self._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata))

    def _text(self, string, point, options):
        self.__objects.append(GraphicPrimitive_Text(string, point, options))
        self._extend_x_axis(2*point[0])
        self._extend_y_axis(2*point[1])

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

    def _add_xy_axes(self, subplot, xmin, xmax, ymin, ymax):
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
	for x in xtslmajor:
            if x == y_axis_xpos:
                continue
	    subplot.text(x, xlabel, format(x), fontsize=5,
                         horizontalalignment="center", verticalalignment="top")

	    subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xltheight],
					color='k',linewidth=0.6))

        for x in xtslminor:
	    subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xstheight],
					color='k', linewidth=0.6))

	#the y-axis ticks and labels
        for y in ytslmajor:
            if y == x_axis_ypos:
                continue
	    subplot.text(ylabel, y, format(y), fontsize=5, verticalalignment="center",
					horizontalalignment="right")

	    subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + yltheight], [y, y],
					color='k', linewidth=0.6))

        for y in ytslminor:
	    subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + ystheight], [y, y],
					color='k',linewidth=0.6))

    def show(self, xmin=None, xmax=None, ymin=None, ymax=None,figsize=DEFAULT_FIGSIZE, filename=None):
        """
	Show a graphics image with default image viewer.
	(Current implementation is hackish)

	EXAMPLES:
	    sage: c = circle((1,1), 1, rgbcolor=(1,0,0))
	    sage: c.save("circ1.png", xmin=-1, xmax=3, ymin=-1, ymax=3)

	    To correct the apect ratio of certain graphics, it is necessary
	    to show with a 'figsize' of square dimensions.

	    sage: c.save("circ1.png", figsize=[5,5], xmin=-1, xmax=3, ymin=-1, ymax=3)

        """
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self.save(filename, xmin, xmax, ymin, ymax, figsize)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(sage.misc.log.browser(), filename))
        #os.system('gqview %s >/dev/null&'%filename)

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

    def save(self, filename, xmin=None, xmax=None, ymin=None, ymax=None, figsize=DEFAULT_FIGSIZE,
		fig=None, sub=None, savenow=True):
        """
	Save the graphics to an image file of type: PNG, PS, or SVG,
	depending on the file extension you give the filename.
        Extension types can be: '.png', '.ps', '.svg'

	EXAMPLES:
	    sage: c = circle((1,1),1,rgbcolor=(1,0,0))
	    sage: c.save("circle.png", xmin=-1,xmax=3,ymin=-1,ymax=3)

	    To correct the apect ratio of certain graphics, it is necessary
	    to show with a 'figsize' of square dimensions.

	    sage: c.save("circle.png", figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)

	"""

        xmin,xmax,ymin,ymax = self._prepare_axes(xmin, xmax, ymin, ymax)

	figure = fig
	if not figure:
             figure = Figure(figsize)
	subplot = sub
	if not subplot:
            subplot = figure.add_subplot(111)
        subplot.xaxis.set_visible(False)
       	subplot.yaxis.set_visible(False)
        subplot._frameon = False
        subplot.set_xlim(xmin, xmax)
        subplot.set_ylim(ymin, ymax)
       	self._add_xy_axes(subplot, xmin, xmax, ymin, ymax)

	for g in self.__objects:
	    g._render_on_subplot(subplot)

        # you can output in PNG, PS, or SVG format, depending on the file extension
	if savenow:
	    try:
	        ext = filename.split('.')[1].lower()
	    except IndexError:
		print "file type must be either 'png' or 'ps' or 'svg'"
	    if ext == 'ps':
                canvas = FigureCanvasPS(figure)
	    elif ext == 'svg':
                canvas = FigureCanvasSVG(figure)
	    elif ext == 'png':
                canvas = FigureCanvasAgg(figure)
	    else:
	        raise ValueError, "file type must be either 'png' or 'ps' or 'svg'"
            canvas.print_figure(filename, dpi=80)

################## Graphics Primitives ################

class GraphicPrimitive:
    def __repr__(self):
        return "Graphics primitive"


class GraphicPrimitive_Line(GraphicPrimitive):
    """
    Primitive class that initializes the
    line graphics type
    """
    def __init__(self, xdata, ydata, options):
        self.xdata = xdata
        self.ydata = ydata
        self.options = options

    def __repr__(self):
        return "Line defined by %s points"%len(self)

    def __getitem__(self, i):
        return self.xdata[int(i)], self.ydata[int(i)]

    def __setitem__(self, i, point):
        i = int(i)
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        return len(self.xdata)

    def __append__(self, point):
        self.xdata.append(float(point[0]))
        self.ydata.append(float(point[1]))

    def _render_on_subplot(self, subplot):
        options = self.options
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
        self.options = options

    def __repr__(self):
        return "Circle defined by (%s,%s) with r=%s"%(self.x, self.y, self.r)

    def _render_on_subplot(self, subplot):
        options = self.options
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

class GraphicPrimitive_Disk(GraphicPrimitive):
    """
    Primitive class that initializes the
    disk graphics type
    """
    def __init__(self, point, r, theta1, theta2, options):
        self.x = point[0]
        self.y = point[1]
	self.r = r
	self.theta1 = theta1
	self.theta2 = theta2
	self.options = options

    def __repr__(self):
        return "Disk defined by (%s,%s) with r=%s with theta (%s, %s)"%(self.x,
			 self.y, self.r, self.theta1, self.theta2)

    def _render_on_subplot(self, subplot):
        options = self.options
        p = patches.Wedge((float(self.x), float(self.y)), float(self.r), float(self.theta1),
 				float(self.theta2))
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
        self.options = options

    def __repr__(self):
        return "Point set defined by %s point(s)"%len(self.xdata)

    def __getitem__(self, i):
	if i == 0:
            return self.xdata
	elif i == 1:
	    return self.ydata
	else:
	    return IndexError("Index out of range")

    def __setitem__(self, i, val):
	i = int(i)
       	if i == 0:
            self.xdata = float(val)
	elif i == 1:
	    self.ydata = float(val)
	else:
	    return IndexError("Index out of range")

    def _render_on_subplot(self,subplot):
        options = self.options
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
        self.options = options

    def __repr__(self):
        return "Polygon defined by %s points"%len(self)

    def __getitem__(self, i):
        return self.xdata[i], self.ydata[i]

    def __setitem__(self, i, point):
        i = int(i)
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        return len(self.xdata)

    def __append__(self, point):
        self.xdata.append(float(point[0]))
        self.ydata.append(float(point[1]))

    def _render_on_subplot(self, subplot):
        options = self.options
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
    Primitive class that initializes the
    text graphics type

    """
    def __init__(self, string, point, options):
	self.string = string
        self.x = point[0]
        self.y = point[1]
        self.options = options

    def __repr__(self):
        return "%s at the point (%s,%s)"%(self.string, self.x, self.y)

    def _render_on_subplot(self, subplot):
        options = self.options
	c = to_mpl_color(options['rgbcolor'])
	f = int(options['fontsize'])
        va = options['vertical_alignment']
	ha = options['horizontal_alignment']
	subplot.text(float(self.x), float(self.y), self.string, color=c, fontsize=f,
					verticalalignment=va,horizontalalignment=ha)


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

class GraphicPrimitiveFactory_disk(GraphicPrimitiveFactory):
    def __call__(self, point, radius, theta1, theta2, **kwds):
	options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata((float(point[0]), float(point[1])),float(radius),
				float(theta1), float(theta2), options=options)

class GraphicPrimitiveFactory_text(GraphicPrimitiveFactory):
    def __call__(self, string, point, **kwds):
	options = dict(self.options)
	for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata(string, (float(point[0]), float(point[1])), options=options)

class GraphicPrimitiveFactory_point(GraphicPrimitiveFactory):
    def __call__(self, point, **kwds):
	options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        return self._from_xdata_ydata(float(point[0]), float(point[1]), options=options)

class GraphicPrimitiveFactory_from_point_list(GraphicPrimitiveFactory):
    def __call__(self, points, coerce=True, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
	if not isinstance(points[0], (list, tuple)):
	    xdata = [float(points[0])]
	    ydata = [float(points[1])]
	else:
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
    ({\em much} better that the gnuplot version!):

        sage: L = [[(pari(1/2 + i*I/10).zeta().real()).precision(1),(pari(1/2 + i*I/10).zeta().imag()).precision(1)] for i in range (0,300)]
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
        pass

    def __repr__(self):
        return "type line? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g._line(xdata, ydata, options)
        return g

# unique line instance
line = LineFactory()
line.options['alpha'] = 1
line.options['rgbcolor'] = (0,0,0)
line.options['thickness'] = 1


class CircleFactory(GraphicPrimitiveFactory_circle):
    """
    A circle at a point = (x,y) with radius = r
    Type circle.options to see all options

    EXAMPLES:
    sage: c = circle((1,1),1,rgbcolor=(1,0,0))
    sage: c.save("circ1.png", xmin=-1,xmax=3,ymin=-1,ymax=3)

    To correct the apect ratio of certain graphics, it is necessary
    to show with a 'figsize' of square dimensions.

    sage: c.save("circ1.png", figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)

    Here we make an more complicated plot with many circles of different colors

    sage: g = Graphics()
    sage: step=6; ocur=1/5; paths=16;
    sage: for r in range(1,paths+1):
    ...       for x,y in [((r+ocur)*cos(n), (r+ocur)*sin(n)) for n in srange(0, 2*pi+pi/step, pi/step)]:
    ...           g += circle((x,y), ocur, rgbcolor=hue(r/paths))
    ...       rnext = (r+1)^2
    ...       ocur = (rnext-r)-ocur
    ...
    sage: g.save("c2.png", xmin=-(paths+1)^2, xmax=(paths+1)^2, ymin=-(paths+1)^2, ymax=(paths+1)^2, figsize=[6,6])
    sage: os.unlink('c2.png')     # remove temp file
    """
    def _reset(self):
	self.options={'fill':False,'thickness':1,'rgbcolor':(0, 0, 0),'resolution':40}

    def __repr__(self):
        return "type circle? for help and examples"

    def _from_xdata_ydata(self, point, r, options):
        g = Graphics()
        g._circle(float(point[0]), float(point[1]), float(r), options)
        return g

#an unique circle instance
circle = CircleFactory()
circle.options['alpha'] = 1
circle.options['fill'] = False
circle.options['resolution'] = 40
circle.options['thickness'] = 1

class DiskFactory(GraphicPrimitiveFactory_disk):
    """
    A disk at a point = (x,y) with radius = r from (theta1, theta2)
    Type disk.options to see all options

    EXAMPLES:
    Make some dangerous disks:

	sage: bl = disk((0.0,0.0),1,180,270,rgbcolor=(1,1,0))
	sage: tr = disk((0.0,0.0),1,0,90,rgbcolor=(1,1,0))
	sage: tl = disk((0.0,0.0),1,90,180,rgbcolor=(0,0,0))
	sage: br = disk((0.0,0.0),1,270,360,rgbcolor=(0,0,0))
	sage: P  = tl+tr+bl+br
	sage: P.save("danger.png",figsize=(4,4),xmin=-2,xmax=2,ymin=-2,ymax=2)
        sage: os.unlink('danger.png')     # remove temp file

    """
    def _reset(self):
	self.options={'fill':True,'resolution':40,'thickness':0}

    def __repr__(self):
        return "type disk? for help and examples"

    def _from_xdata_ydata(self, point, r, theta1, theta2, options):
        g = Graphics()
        g._disk((float(point[0]), float(point[1])), float(r), float(theta1), float(theta2), options)
        return g

#an unique disk instance
disk = DiskFactory()
disk.options['alpha'] = 1
disk.options['fill'] = True
disk.options['rgbcolor'] = (0,0,0)
disk.options['thickness'] = 0

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
        self.options = {'pointsize':10,'faceted':False,'rgbcolor':(0,0,0)}

    def __repr__(self):
        return "type point? for options help"

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g._point(xdata, ydata, options)
        return g

# unique point instance
point=PointFactory()
point.options['alpha'] = 1
point.options['pointsize'] = 10
point.options['faceted'] = False

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
        self.options = {'fontsize':10,'rgbcolor':(0,0,0)}

    def __repr__(self):
        return "type text? for help and examples"

    def _from_xdata_ydata(self, string, point, options):
        g = Graphics()
        g._text(string, point, options)
        return g

# unique text instance
text = TextFactory()
text.options['fontsize'] = 10
text.options['horizontal_alignment']= "center"
text.options['vertical_alignment'] = "center"

class PolygonFactory(GraphicPrimitiveFactory_from_point_list):
    """
    Type polygon.options for a dictionary of the default
    options for polygons.  You can change this to change
    the defaults for all future polygons.  Use polygon.reset()
    to reset to the default options.

    EXAMPLES:
    We create a purple-ish polygon:
	sage: polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))
	Graphics object consisting of 1 graphics primitives:
        	0 -- Polygon defined by 3 points

    Some modern art -- a random polygon:
	sage: v = [(randrange(-5,5), randrange(-5,5)) for _ in range(10)]
	sage: polygon(v)
	Graphics object consisting of 1 graphics primitives:
        	0 -- Polygon defined by 10 points

    A purple hexagon:

        sage: L = [[cos(pi*i/3),sin(pi*i/3)] for i in range(6)]
        sage: p = polygon(L, rgbcolor=(1,0,1))

    A green limacon of Pascal:

        sage: L = [[1+cos(pi*i/100)-cos(pi*i/50),sin(pi*i/100)-sin(pi*i/50)] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/4,1,1/2))

    A green deltoid:

        sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))

    A blue figure 8:

        sage: L = [[2*cos(pi*i/100)*sqrt(1-sin(pi*i/100)^2),2*sin(pi*i/100)*sqrt(1-sin(pi*i/100)^2)] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))

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

    A blue 8-leaved petal:

        sage: L = [[sin(5*pi*i/100)^2*cos(pi*i/100)^3,sin(5*pi*i/100)^2*sin(pi*i/100)] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1/3,1/2,3/5))

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
        pass

    def __repr__(self):
        return "SAGE polygon; type polygon? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g._polygon(xdata, ydata, options)
        return g


# unique polygon instance
polygon = PolygonFactory()
polygon.options['alpha'] = 1
polygon.options['rgbcolor'] = (0,0,0)
polygon.options['thickness'] = 0

class PlotFactory(GraphicPrimitiveFactory):
    """
    Type plot.options for a dictionary of the default
    options for plots.  You can change this to change
    the defaults for all future plots.  Use plot.reset()
    to reset to the default options.

    EXAMPLES:
    We draw the graph of an elliptic curve as the union
    of graphs of 2 functions.
        sage: def h1(x): return sqrt(x^3  - 1)
        sage: def h2(x): return -sqrt(x^3  - 1)
	sage: plot([h1, h2], 1,4)
	Graphics object consisting of 2 graphics primitives:
        	0 -- Line defined by 204 points
        	1 -- Line defined by 204 points
    """
    def _reset(self):
        o = self.options
        o['plot_points'] = 200
        o['plot_division'] = 1000      # change when fix adapter
        o['max_bend'] = 0.1            # change when fix adapter

    def __repr__(self):
        return "plot; type plot? for help and examples."

    def __call__(self, f, xmin, xmax, parametric=False, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        if isinstance(f, (list, tuple)) and not parametric:
            G = Graphics()
            for i in range(0, len(f)):
                G += plot(f[i], xmin, xmax, **kwds)
            return G
	if parametric:
	    g,f = f #f=(f,g)
        xmin = float(xmin)
        xmax = float(xmax)
        plot_points = int(options['plot_points'])
        delta = (xmax - xmin) / plot_points
        data = []
        for i in xrange(plot_points + 1):
            x = xmin + i*delta
            data.append((x, float(f(x))))
        # adaptive refinement
        i, j = 0, 0
        max_bend = float(options['max_bend'])
        plot_division = int(options['plot_division'])
        while i < len(data) - 1:
            if abs(data[i+1][1] - data[i][1]) > max_bend:
                x = (data[i+1][0] + data[i][0])/2
                y = float(f(x))
                data.insert(i+1, (x,y))
                j += 1
                if j > plot_division:
                    # wrong -- just for testing (so no infinite loop)
                    break
            else:
                i += 1
	if parametric:
	    data = [(float(g(x)),y) for x,y in data]
        g = line(data, coerce=False, **options)
        return g

# unique plot instance
plot = PlotFactory()

########## misc functions ###################

#parametric plot
def parametric_plot((f,g), xmin, xmax, **kwargs):
    """
    parametric plot takes two functions as a
    list or a tuple and make a plot with the
    first function as the x-value and the
    second function as the y-value

    EXAMPLE:
	f = lambda x: sin(x)
	g = lambda x: sin(2*x)
	p3 = parametric_plot((f,g),0,2*pi,rgbcolor=hue(0.6))

    """

    return plot((f,g), xmin, xmax, parametric=True, **kwargs)

#list plot
def list_plot(data, plotjoined=False, **kwargs):
    """
        list_plot takes a single list of data,
    in which case it forms a list of tuples (i,di)
    where i goes from 0 to len(data)-1 and di is
    the ith data value, and puts points at those
    tuple values
	list_plot also takes a list of tuples
    (dxi, dyi) where dxi is the ith data representing
    the x-value, and dyi is the ith y-value
        if plotjoined=True, then a line spanning
    all the data is drawn instead

    EXAMPLES:
	sage: l = list_plot([i^2 for i in range(5)]); l
	Graphics object consisting of 1 graphics primitives:
        	0 -- Point set defined by 5 point(s)

	sage: r = [(random(),random()) for _ in range(20)]

    Here are a bunch of random red points:
	sage: list_plot(r,rgbcolor=(1,0,0)).save('random.png')
        sage: os.unlink('random.png')

    This gives all the random points joined in a purple line:
	sage: list_plot(r, plotjoined=True, rgbcolor=(1,0,1)).save("randomjoin.png")
    """
    if not isinstance(data[0], (list, tuple)):
	 data = zip(range(len(data)),data)
    if plotjoined:
	return line(data, **kwargs)
    return point(data, **kwargs)

def to_float_list(v):
    return [float(x) for x in v]

def to_mpl_color(c):
    return (float(c[0]), float(c[1]), float(c[2]))

def hue(h, s=1, v=1):
    """
      hue(h,s=1,v=1) where 'h' stands for hue,
      's' stands for saturation, 'v' stands for value.
      hue returns a list of rgb intensities (r, g, b)
      All values are in range 0 to 1.

      EXAMPLES:
	sage: hue(0.6)
        (0.0, 0.40000000000000036, 1.0)

	hue is an easy way of getting a broader
	range of colors for graphics

	sage: p = plot(sin, -2, 2, rgbcolor=hue(0.6))

    """
    if h > 1:
        h = modf(h)[0]
    if s > 1:
        s = modf(s)[0]
    if v > 1:
        v = modf(v)[0]
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

    def __append__(self, g):
        self._glist.append(g)

    def _render(self, filename):
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
	figure = Figure()
	for i,g in zip(range(1, dims+1), glist):
            subplot = figure.add_subplot(rows, cols, i)
	    g.save(filename, fig=figure, sub=subplot, savenow = (i==dims))   # only save if i==dims.

    def save(self, filename="ga.png"):
	"""
	save the \code{graphics_array} to
        (for now) a png called 'filename'.
	"""
	self._render(filename)

    def show(self, filename=None):
	"""
	show the \code{graphics_array} in
        the users browser.
	"""
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self.render(filename)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(
                         sage.misc.log.browser(), filename))
        #os.system('gqview %s >/dev/null&'%filename)

def graphics_array(array,filename="ga.png"):
    r"""
    \code{graphics_array} take a list of lists (or tuples)
    of graphics objects and plots them all on one canvas (single plot).

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
    G = GraphicsArray(array)
    return G






