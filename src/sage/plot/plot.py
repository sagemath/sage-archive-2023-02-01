"""
Plotting

AUTHORS:
    -- Alex Clemesha and William Stein (2006-04-10): initial version
    -- David Joyner: examples

EXAMPLES:
    sage: G = plot(cos, -5, 5, thickness=5, rgbcolor=(0.5,1,0.5))
    sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))

Construct reflection of above polygon about the y-axis by iterating
over list of first-coordinates of the first graphic element of P
(which is the actual Polygon -- note that P is a Graphics object,
which consists of a single polygon):

    sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
    sage: H = G + P + Q
    sage: H
    Graphics object consiting of 3 graphics primitives:
            0 -- Line defined by 78 points
            1 -- Polygon defined by 3 points
            2 -- Polygon defined by 3 points
    sage: type(H)
    <class 'sage.plot.plot.Graphics'>
    sage: H[1]
    Polygon defined by 3 points
    sage: list(H[1])
    [(1.0, 2.0), (5.0, 6.0), (5.0, 0.0)]
"""

#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

from sage.ext.sage_object import SageObject

try:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.transforms import Value
    from matplotlib.figure import Figure
    import matplotlib.patches as patches
except ImportError:
    pass

from axes import find_axes

class Graphics(SageObject):
    def __init__(self):
        self.__xmin = -1
        self.__xmax = 1
        self.__ymin = -1
        self.__ymax = 1
        self.__objects = []

    def xmax(self):
        return self.__xmax

    def xmin(self):
        return self.__xmin

    def ymax(self):
        return self.__ymax

    def ymin(self):
        return self.__ymin

    def _repr_(self):
        pr = ''
        i = 0
        for x in self:
            pr += '\n\t%s -- %s'%(i, x)
            i += 1
        return "Graphics object consiting of %s graphics primitives:%s"%(
            len(self), pr)

    def __getitem__(self, i):
        return self.__objects[int(i)]

    def __len__(self):
        return len(self.__objects)

    def __delitem__(self, i):
        del self.__objects[int(i)]

    def __setitem__(self, i, x):
        if not isinstance(x, GraphicPrimitive):
            raise TypeError, "x must be a GraphicPrimitive"
        self.__objects[int(i)] = x


    def __add__(self, other):
        """
        EXAMPLE:
            sage: h1 = lambda x : sqrt(x^3  - 1)
            sage: h2 = lambda x : -sqrt(x^3  - 1)
            sage: g1 = plot(h1, 1, 5)
            sage: g2 = plot(h2, 1, 5)
            sage: h = g1 + g2; h
            Graphics object consiting of 2 graphics primitives:
                0 -- Line defined by 161 points
                1 -- Line defined by 161 points
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

    def line(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Line(xdata, ydata, options))
        self._extend_x_axis(min(xdata))
        self._extend_x_axis(max(xdata))
        self._extend_y_axis(min(ydata))
        self._extend_y_axis(max(ydata))

    def polygon(self, xdata, ydata, options):
        self.__objects.append(GraphicPrimitive_Polygon(xdata, ydata, options))
        self._extend_x_axis(min(xdata))
        self._extend_x_axis(max(xdata))
        self._extend_y_axis(min(ydata))
        self._extend_y_axis(max(ydata))

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

    def _add_xy_axes(self, subplot, xmin, xmax, ymin, ymax):
        yspan = ymax - ymin
        xspan = xmax - xmin

	#evalute find_axes for x values and y ticks
	y_axis_xpos, xstep, xtslminor, xtslmajor = find_axes(xmin, xmax)
	yltheight = 0.02 * xspan
	ystheight = 0.4  * yltheight
	ylabel    = y_axis_xpos -2*ystheight

	#evalute find_axes for y values and x ticks
	x_axis_ypos, ystep, ytslminor, ytslmajor = find_axes(ymin, ymax)
	xltheight = 0.02 * yspan
        xstheight = 0.4  * xltheight
	xlabel    = x_axis_ypos - xltheight

	#the x axis line
        subplot.add_line(patches.Line2D([xmin, xmax], [x_axis_ypos, x_axis_ypos],
                                        color='k', linewidth=0.8))

	#the y axis line
        subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos],[ymin, ymax],
                                        color='k', linewidth=0.8))


        def format(z):
            s = str(z)
            if s[-2:] == '.0':
                return s[:-2]
            return s


	#the x-axis ticks and labels
	for x in xtslmajor:
            if x == y_axis_xpos:
                continue
	    subplot.text(x, xlabel, format(x), fontsize=6,
                         horizontalalignment="center", verticalalignment="top")

	    subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xltheight],
					color='k',linewidth=0.8))

        for x in xtslminor:
	    subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xstheight],
					color='k', linewidth=0.8))

	#the y-axis ticks and labels
        for y in ytslmajor:
            if y == x_axis_ypos:
                continue
	    subplot.text(ylabel, y, format(y), fontsize=6, verticalalignment="center",
					horizontalalignment="right")

	    subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + yltheight], [y, y],
					color='k', linewidth=0.8))

        for y in ytslminor:
	    subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + ystheight], [y, y],
					color='k',linewidth=0.8))

    def show(self, xmin=None, xmax=None, ymin=None, ymax=None):
        # Todo -- mainly for testing.
        t = '/tmp/sage.png'
        self.save(t, xmin, xmax, ymin, ymax)
        os.system('gqview %s >/dev/null&'%t)

    def save(self, filename, xmin=None, xmax=None,
             ymin=None, ymax=None):

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

        figure = Figure()
        subplot = figure.add_subplot(111)
        subplot.xaxis.set_visible(False)
        subplot.yaxis.set_visible(False)
        subplot._frameon = False
        subplot.set_xlim(xmin, xmax)
        subplot.set_ylim(ymin, ymax)

        self._add_xy_axes(subplot, xmin, xmax, ymin, ymax)

        for g in self.__objects:
            p = g.patch()
            if isinstance(p, patches.Polygon):
                subplot.add_patch(p)
            else:
                subplot.add_line(p)

        # this is the only reason it outputs png...
        canvas = FigureCanvasAgg(figure)
        canvas.print_figure(filename)


################## Graphics Primitives ################

class GraphicPrimitive:
    def __repr__(self):
        return "Graphics primitive"


class GraphicPrimitive_Line(GraphicPrimitive):
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

    def patch(self):
        options = self.options
        p = patches.Line2D(self.xdata, self.ydata)
        p.set_linewidth(float(options['thickness']))
        p.set_color(to_mpl_color(options['rgbcolor']))
        return p


class GraphicPrimitive_Polygon(GraphicPrimitive):
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

    def patch(self):
        options = self.options
        p = patches.Polygon([(self.xdata[i],self.ydata[i]) for i in xrange(len(self.xdata))])
        p.set_linewidth(float(options['thickness']))
        c = to_mpl_color(options['rgbcolor'])
        # I think that in Mathematica there is no way to
        # color just the edges of the polygon...
        p.set_edgecolor(c)
        p.set_facecolor(c)
        return p

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
        self.options = {'thickness':1,
                        'rgbcolor':(0,0,0)}
        self._reset()

    def _coerce(self, xdata, ydata):
        return to_float_list(xdata), to_float_list(ydata)


class GraphicPrimitiveFactory_from_point_list(GraphicPrimitiveFactory):
    def __call__(self, points, coerce=True, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v

        if coerce:
            xdata = []
            ydata = []
            for z in points:
                xdata.append(float(z[0]))
                ydata.append(float(z[1]))
        else:
            xdata = [z[0] for z in points]
            ydata = [z[1] for z in points]

        return self._from_xdata_ydata(xdata, ydata, coerce=False, options=options)

class Line(GraphicPrimitiveFactory_from_point_list):
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
        return "SAGE line plotter; type line? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g.line(xdata, ydata, options)
        return g

# unique plot instance
line = Line()


class Polygon(GraphicPrimitiveFactory_from_point_list):
    """
    Type polygon.options for a dictionary of the default
    options for polygons.  You can change this to change
    the defaults for all future polygons.  Use polygon.reset()
    to reset to the default options.

    EXAMPLES:
    We create a purple-ish polygon:
        sage: polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))
        Graphics object consiting of 1 graphics primitives:
            0 -- Polygon defined by 3 points


    Some modern art -- a random polygon:
        sage: v = [(randrange(-5,5), randrange(-5,5)) for _ in range(10)]
        sage: polygon(v)
        Graphics object consiting of 1 graphics primitives:
            0 -- Polygon defined by 10 points

    Here's a yellow circle:

        sage: L = [[cos(pi*i/100),sin(pi*i/100)] for i in range(200)]
        sage: p = polygon(L, rgbcolor=(1,1,0))

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
        -- David Joyner (2060-04-14): the long list of examples above.
    """
    def _reset(self):
        pass

    def __repr__(self):
        return "SAGE polygon plotter; type polygon? for help and examples."

    def _from_xdata_ydata(self, xdata, ydata, coerce, options):
        if coerce:
            xdata, ydata = self._coerce(xdata, ydata)
        g = Graphics()
        g.polygon(xdata, ydata, options)
        return g


# unique plot instance
polygon = Polygon()


class Plot(GraphicPrimitiveFactory):
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
        Graphics object consiting of 2 graphics primitives:
            0 -- Line defined by 103 points
            1 -- Line defined by 103 points
    """
    def _reset(self):
        o = self.options
        o['plot_points'] = 25
        o['plot_division'] = 1000      # change when fix adapter
        o['max_bend'] = 0.1            # change when fix adapter

    def __repr__(self):
        return "SAGE Plotter; type plot? for help and examples."

    def __call__(self, f, xmin, xmax, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        if isinstance(f, (list, tuple)):
            i = 0
            G = plot(f[0], xmin, xmax, **kwds)
            for i in range(1, len(f)):
                G += plot(f[i], xmin, xmax, **kwds)
            return G

        xmin = float(xmin)
        xmax = float(xmax)
        plot_points = int(options['plot_points'])
        delta = (xmax - xmin) / plot_points
        data = []
        for i in xrange(plot_points):
            x = xmin + i*delta
            data.append((x, float(f(x))))

        # adaptive refinement
        i = 0
        j = 0
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

        g = line(data, coerce=False, **options)

        return g

# unique plot instance
plot = Plot()




########## misc functions ###################

def to_float_list(v):
    return [float(x) for x in v]

def to_mpl_color(c):
    return (float(c[0]), float(c[1]), float(c[2]))
