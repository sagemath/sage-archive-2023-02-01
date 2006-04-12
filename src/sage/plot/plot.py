"""
Plotting

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

import os

from sage.ext.sage_object import SageObject

try:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.transforms import Value
    from matplotlib.figure import Figure
    import matplotlib.patches as patches
except ImportError:
    pass

class Graphics(SageObject):
    def __init__(self):
        self.__xmin = None
        self.__xmax = None
        self.__ymin = None
        self.__ymax = None
        self.__objects = []

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
        if 0 < xmin:
            xat = xmin  + 0.05 * (xmax - xmin)
        else:
            xat = 0
        if 0 < ymin:
            yat = ymin  + 0.05 * (ymax - ymin)
        else:
            yat = 0
        subplot.add_line(patches.Line2D([xmin, xmax], [yat,yat],
                                        color='k', linewidth=1))
        subplot.add_line(patches.Line2D([xat,xat], [ymin, ymax],
                                        color='k', linewidth=1))
        # todo: add tick marks and text like in mathematica

    def show(self, xmin=None, xmax=None, ymin=None, ymax=None):
        # Todo -- mainly for testing.
        t = '/tmp/sage.png'
        self.save(t, xmin, xmax, ymin, ymax)
        os.system('gqview %s &'%t)

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
    """
    Type line.options for a dictionary of the default
    options for lines.  You can change this to change
    the defaults for all future lines.  Use line.reset()
    to reset to the default options.
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
        o['plot_division'] = 1000000   # change when fix adapter
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
