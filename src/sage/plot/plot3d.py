
import os


import sage.misc.viewer
import sage.misc.misc
from sage.ext.sage_object import SageObject

from mpl3d import mplot3d

from plot import to_float_list, to_mpl_color
import plot

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_ps import FigureCanvasPS
from matplotlib.backends.backend_svg import FigureCanvasSVG

import matplotlib.numerix as nx

class Graphics3d(SageObject):
    def __init__(self):
        #print "WARNING: 3d Graphics -- WORK IN PROGRESS!!"
        self.__xmin = -1
        self.__xmax = 1
        self.__ymin = -1
        self.__ymax = 1
        self.__zmin = -1
        self.__zmax = 1
        self.xlabel = 'x'
        self.ylabel = 'y'
        self.zlabel = 'z'
        self.__objects = []

    def xmax(self):
        """
	sage: G = Graphics3d(); G
	3d Graphics object consisting of 0 graphics primitives:
	sage: G.xmax()
	1
        """
        return self.__xmax

    def xmin(self):
	"""
	sage: G = Graphics3d(); G
	3d Graphics object consisting of 0 graphics primitives:
	sage: G.xmin()
	-1
	"""
        return self.__xmin

    def ymax(self):
	"""
	sage: G = Graphics3d(); G
	3d Graphics object consisting of 0 graphics primitives:
	sage: G.ymax()
	1
	"""
        return self.__ymax

    def ymin(self):
	"""
	sage: G = Graphics3d(); G
	3d Graphics object consisting of 0 graphics primitives:
	sage: G.ymin()
	-1
	"""
        return self.__ymin

    def zmax(self):
	"""
	sage: G = Graphics3d(); G
	3d Graphics object consisting of 0 graphics primitives:
	sage: G.ymax()
	1
	"""
        return self.__zmax

    def zmin(self):
	"""
	sage: G = Graphics3d(); G
        3d Graphics object consisting of 0 graphics primitives:
	sage: G.ymin()
	-1
	"""
        return self.__zmin

    def _repr_(self):
        pr, i = '', 0
        for x in self:
            pr += '\n\t%s -- %s'%(i, x)
            i += 1
        return "3d Graphics object consisting of %s graphics primitives:%s"%(
            len(self), pr)

    def __getitem__(self, i):
	"""
	Returns the i-th graphics primitive object:
	"""
        return self.__objects[int(i)]

    def __len__(self):
        """
        If G is of type Graphics3d, then len(G) gives the number of
	distinct graphics primitives making up that object.

	EXAMPLES:
        """
        return len(self.__objects)

    def __delitem__(self, i):
        """
        If G is of type Graphics, then del(G[i]) removes the ith
	distinct graphic primitive making up that object.
        """
        del self.__objects[int(i)]

    def __setitem__(self, i, x):
        """
	You can replace an GraphicsPrimitive3d (point, line, circle, etc...)
	in a Graphics3d object G with any other GraphicsPrimitive3d.
        """
        if not isinstance(x, Graphic3dPrimitive):
            raise TypeError, "x must be a Graphic3dPrimitive"
        self.__objects[int(i)] = x

    def __radd__(self, other):
        if isinstance(other, int) and other == 0:
            return self
        raise TypeError

    def __add__(self, other):
        """
        If you have any Graphics3d object G1, you can always add any
	other amount of Graphics3d objects G2,G3,...  to form a new
	Graphics object: G4 = G1 + G2 + G3
        """
        if isinstance(other, int) and other == 0:
            return self
        if not isinstance(other, Graphics3d):
            raise TypeError, "other (=%s) must be a Graphics objects"%other
        g = Graphics3d()
        g.__xmin = min(self.__xmin, other.__xmin)
        g.__xmax = max(self.__xmax, other.__xmax)
        g.__ymin = min(self.__ymin, other.__ymin)
        g.__ymax = max(self.__ymax, other.__ymax)
        g.__zmin = min(self.__zmin, other.__zmin)
        g.__zmax = max(self.__zmax, other.__zmax)
        g.__objects = self.__objects + other.__objects
        return g

    def _point3d(self, xdata, ydata, zdata, options):
        self.__objects.append(Graphic3dPrimitive_Point(xdata, ydata, zdata, options))
        try:
            self._extend_axes(min(xdata), max(xdata), min(ydata), max(ydata), min(zdata), max(zdata))
        except ValueError:
            pass

    def _line3d(self, xdata, ydata, zdata, options):
        self.__objects.append(Graphic3dPrimitive_Line(xdata,
                                                      ydata,
                                                      zdata, options))
        try:
            self._extend_axes(min(xdata), max(xdata), min(ydata),
                              max(ydata), min(zdata), max(zdata))
        except ValueError:
            pass

    def _surface3d(self, data, options):
        self.__objects.append(Graphic3dPrimitive_Surface(data, options))
        try:
            xdata = sum([[z[0] for z in w] for w in data], [])
            ydata = sum([[z[1] for z in w] for w in data], [])
            zdata = sum([[z[2] for z in w] for w in data], [])
            self._extend_axes(min(xdata), max(xdata), min(ydata),
                              max(ydata), min(zdata), max(zdata))
        except ValueError:
            pass


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

    def _extend_z_axis(self, z):
        zmin = self.__zmin
        zmax = self.__zmax
        if zmin is None or z < zmin:
            self.__zmin = z
        elif zmax is None or z > zmax:
            self.__zmax = z

    def _extend_axes(self, xmin, xmax, ymin, ymax, zmin, zmax):
    	self._extend_x_axis(xmin)
    	self._extend_x_axis(xmax)
    	self._extend_y_axis(ymin)
    	self._extend_y_axis(ymax)
    	self._extend_z_axis(zmin)
    	self._extend_z_axis(zmax)

    def show(self,
             xmin=None, xmax=None,
             ymin=None, ymax=None,
             zmin=None, zmax=None,
             figsize  = plot.DEFAULT_FIGSIZE,
             filename = None,
             dpi      = None):
        """
	Show this Graphics3d image with the default image viewer.
        """
        if plot.EMBEDDED_MODE:
            self.save(filename, xmin, xmax, ymin, ymax, zmin, zmax, figsize, dpi=dpi)
            return
        if filename is None:
            filename = sage.misc.misc.tmp_filename() + '.png'
        self.save(filename, xmin, xmax, ymin, ymax, zmin, zmax, figsize, dpi=dpi)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(sage.misc.viewer.browser(), filename))

    def _prepare_axes(self, xmin, xmax, ymin, ymax, zmin, zmax):
        if self.__xmin is None:
            self.__xmin, self.__xmax, self.__ymin, self.__ymax, self.__zmin, self.__zmax = 0, 0, 0, 0, 0, 0

        if xmin is None:
            xmin = self.__xmin
        if xmax is None:
            xmax = self.__xmax
        if ymin is None:
            ymin = self.__ymin
        if ymax is None:
            ymax = self.__ymax
        if zmin is None:
            zmin = self.__zmin
        if zmax is None:
            zmax = self.__zmax

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

        if zmax < zmin:
            zmax, zmin = zmin, zmax
        elif zmax == zmin:
            z = zmax
            if z == 0:
                z = 1
            zmax = 2*z
            zmin = 0

        xmin -= 0.1*(xmax-xmin)
        xmax += 0.1*(xmax-xmin)
        ymin -= 0.1*(ymax-ymin)
        ymax += 0.1*(ymax-ymin)
        zmin -= 0.1*(zmax-zmin)
        zmax += 0.1*(zmax-zmin)

	return xmin,xmax,ymin,ymax,zmin,zmax

    def save(self, filename=None,
             xmin=None, xmax=None,
             ymin=None, ymax=None,
             zmin=None, zmax=None,
             figsize=plot.DEFAULT_FIGSIZE,
             fig=None, sub=None,
             savenow=True, dpi=None,
             elevation = 30, azimuth=-60):
        """
	Save this Graphics3d object to an image file of type: png, eps,
	or svg, depending on the file extension you give the filename.
	Extension types can be: '.png', '.eps', '.svg'.
	"""
        if filename is None:
            i = 0
            while os.path.exists('sage%s.png'%i):
                i += 1
            filename = 'sage%s.png'%i

        try:
            ext = os.path.splitext(filename)[1].lower()
        except IndexError:
            raise ValueError, "file extension must be either 'png', 'ps', 'eps', 'svg' or 'sobj'"

        if ext == '' or ext == '.sobj':
            SageObject.save(self, filename)
            return

        xmin,xmax,ymin,ymax,zmin,zmax = self._prepare_axes(xmin, xmax, ymin, ymax, zmin, zmax)

	figure = fig
	if not figure:
            figure = Figure(figsize)

	subplot = sub
	if not subplot:
            # we call it a subplot for consistency with plot 2d.
            # it's really just an axis.
            subplot = mplot3d.Axes3D(figure,
                                     elev = elevation,
                                     azim = azimuth)

        subplot.set_w_xlim(xmin, xmax)
        subplot.set_w_ylim(ymin, ymax)
        subplot.set_w_zlim(zmin, zmax)
        subplot.set_xlabel(str(self.xlabel))
        subplot.set_ylabel(str(self.ylabel))
        subplot.set_zlabel(str(self.zlabel))

	for g in self.__objects:
	    g._render_on_subplot(subplot)

        figure.add_axes(subplot)

        # you can output in PNG, PS, or SVG format, depending on the file extension
	if savenow:
            if ext in ['.ps', '.eps']:
                canvas = FigureCanvasPS(figure)
		if dpi is None:
		    dpi = 72
	    elif ext == '.svg':
                canvas = FigureCanvasSVG(figure)
		if dpi is None:
		    dpi = 80
	    elif ext == '.png':
                canvas = FigureCanvasAgg(figure)
		if dpi is None:
		    dpi = 150
            else:
                raise ValueError, "file extension must be either 'png', 'ps', 'svg' or 'sobj'"
            canvas.print_figure(filename, dpi=dpi)



################## Graphics Primitives ################

class Graphic3dPrimitive(SageObject):
    def _repr_(self):
        return "Graphics primitive 3d"

class Graphic3dPrimitive_Point(Graphic3dPrimitive):
    """
    Primitive class that initializes the
    point graphics type
    """
    def __init__(self, xdata, ydata, zdata, options):
        #see top of this file for Point info
        self.xdata = xdata
        self.ydata = ydata
        self.zdata = zdata
        self.options = options

    def _repr_(self):
        return "Point set defined by %s point(s)"%len(self.xdata)

    def __getitem__(self, i):
	if i == 0:
            return self.xdata
	elif i == 1:
	    return self.ydata
	elif i == 2:
	    return self.zdata
	else:
	    return IndexError("Index out of range")

    def _render_on_subplot(self, subplot):
        options = self.options
	c = to_mpl_color(options['rgbcolor'])
	a = float(options['alpha'])
	s = int(options['pointsize'])
	faceted = options['faceted'] #faceted=True colors the edge of point
        #subplot.scatter3D(self.xdata, self.ydata, self.zdata, s, c, alpha=a, faceted=False)
        subplot.scatter3D(self.xdata, self.ydata, self.zdata, c=c)
        from pylab import array
        subplot.plot3D(array([0,1,0.5]), array([0,1,1.1]), array([0,1,-2.1]))


class Graphic3dPrimitive_Line(Graphic3dPrimitive):
    """
    Primitive class that initializes the
    line graphics type
    """
    def __init__(self, xdata, ydata, zdata, options):
        self.xdata = xdata
        self.ydata = ydata
        self.zdata = zdata
        if not len(xdata) == len(ydata) and len(ydata) == len(zdata):
            raise ValueError, "xdata, ydata, and zdata must all have the same length"
        self.options = options

    def _repr_(self):
        return "Line defined by %s points"%len(self)

    def __getitem__(self, i):
        return self.xdata[int(i)], self.ydata[int(i)], self.zdata[int(i)]

    def __setitem__(self, i, point):
        del self._mpl_data
        i = int(i)
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])
        self.zdata[i] = float(point[2])

    def __len__(self):
        return len(self.xdata)

    def __append__(self, point):
        self.xdata.append(float(point[0]))
        self.ydata.append(float(point[1]))
        self.zdata.append(float(point[2]))



    def _render_on_subplot(self, subplot):
        try:
            xs, ys, zs = self._mpl_data
        except AttributeError:
            xs = nx.array([float(x) for x in self.xdata])
            ys = nx.array([float(y) for y in self.ydata])
            zs = nx.array([float(z) for z in self.zdata])
            self._mpl_data = xs, ys, zs
        P = subplot.plot3D(xs, ys, zs)

        options = self.options
	a = float(options['alpha'])
        for p in P:
            p.set_alpha(a)
            p.set_linewidth(float(options['thickness']))
            p.set_color(to_mpl_color(options['rgbcolor']))




class Graphic3dPrimitive_Surface(Graphic3dPrimitive):
    """
    Primitive class that initializes the surface graphics type

    Imagine a surface embedded in 3 space that has a grid drawn on it.
    The (i,j) position on the grid has coordinates data[i][j].
    """
    def __init__(self, data, options):
        self.data = data
        self.options = options

    def _repr_(self):
        return "Surfce defined by a %s x %s mesh"%(len(self.data), len(self.data[0]))

    def __getitem__(self, t):
        i, j = t
        return self.data[int(i)][int(j)]

    def __setitem__(self, t, point):
        i, j = t
        del self._mpl_data
        self.data[int(i)][int(j)] = point

    def _render_on_subplot(self, subplot):
        try:
            X, Y, Z = self._mpl_data
        except AttributeError:
            n = len(self.data)
            m = len(self.data[0])
            X = nx.zeros((n,m))
            Y = nx.zeros((n,m))
            Z = nx.zeros((n,m))
            for i in range(n):
                for j in range(m):
                    P = self.data[i][j]
                    X[i,j] = P[0]
                    Y[i,j] = P[1]
                    Z[i,j] = P[2]
            self._mpl_data = X,Y,Z

        P = subplot.plot_surface(X,Y,Z, div=int(self.options['div']))



######################################################################
#                                                                    #
#    Graphics3d Primitives Factories -- construct GraphicPrimitives  #
#                                                                    #
######################################################################

class Graphic3dPrimitiveFactory:
    def __init__(self):
        # options for this specific graphics primitive.
        self.reset()

    def reset(self):
        # First the default options for all graphics primitives
        self.options = {'alpha':1,'thickness':1,'rgbcolor':(0,0,0)}
        self._reset()

    def _coerce(self, xdata, ydata, zdata):
        return to_float_list(xdata), to_float_list(ydata), to_float_list(zdata)


class Graphic3dPrimitiveFactory_from_point_list(Graphic3dPrimitiveFactory):
    def __call__(self, points, coerce=True, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v
        if len(points) > 0 and not isinstance(points[0], (list, tuple)):
	    xdata = [float(points[0])]
	    ydata = [float(points[1])]
            zdata = [float(points[2])]
	else:
            if coerce:
                xdata = []
                ydata = []
                zdata = []
                for z in points:
                    xdata.append(float(z[0]))
                    ydata.append(float(z[1]))
                    zdata.append(float(z[2]))
            else:
                xdata = [z[0] for z in points]
                ydata = [z[1] for z in points]
                zdata = [z[2] for z in points]

        return self._from_xdata_ydata_zdata(xdata, ydata, zdata, True, options=options)

class Point3dFactory(Graphic3dPrimitiveFactory_from_point_list):
    def _reset(self):
        self.options = {'alpha':1, 'pointsize':10,'faceted':False,'rgbcolor':(0,0,0)}

    def _repr_(self):
        return "type point? for options help"

    def _from_xdata_ydata_zdata(self, xdata, ydata, zdata, coerce, options):
        if coerce:
            xdata, ydata, zdata = self._coerce(xdata, ydata, zdata)
        g = Graphics3d()
        g._point3d(xdata, ydata, zdata, options)
        return g

# unique point instance
point3d = Point3dFactory()


class Line3dFactory(Graphic3dPrimitiveFactory_from_point_list):
    r"""
    Create the line through the given list of points.
    """
    def _reset(self):
        self.options = {'alpha':1, 'rgbcolor':(0,0,0),
                        'thickness':1}

    def _repr_(self):
        return "type line? for help and examples."

    def _from_xdata_ydata_zdata(self, xdata, ydata, zdata, coerce, options):
        if coerce:
            xdata, ydata, zdata = self._coerce(xdata, ydata, zdata)
        g = Graphics3d()
        g._line3d(xdata, ydata, zdata, options)
        return g


# unique line instance
line3d = Line3dFactory()



class SurfaceFactory(Graphic3dPrimitiveFactory):
    def __call__(self, data, coerce=True, **kwds):
        options = dict(self.options)
        for k, v in kwds.iteritems():
            options[k] = v

        if coerce:
            data = [[[float(a[0]),float(a[1]),float(a[2])]
                              for a in r] for r in data]

        g = Graphics3d()
        g._surface3d(data, options)
        return g

    def _reset(self):
        pass

# unique surface factory
surface = SurfaceFactory()
