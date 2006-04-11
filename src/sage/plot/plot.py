"""
Plotting
"""

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
        self.objects = []

    def _repr_(self):
        return "Graphics object"

    def __add__(self, other):
        """
        EXAMPLE:
            sage: h1 = lambda x : sqrt(x^3  - 1)
            sage: h2 = lambda x : -sqrt(x^3  - 1)
            sage: g1 = plot(h1, 1, 5)
            sage: g2 = plot(h2, 1, 5)
            sage: h = g1 + g2; h
            Graphics object
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
        g.objects = self.objects + other.objects
        return g

    def line2d(self, xdata, ydata):
        p = patches.Line2D(xdata, ydata)
        self.objects.append(p)
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
                                        color='k', linewidth=0.5))
        subplot.add_line(patches.Line2D([xat,xat], [ymin, ymax],
                                        color='k', linewidth=0.5))
        # todo: add tick marks and text like in mathematica


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

        for g in self.objects:
            subplot.add_line(g)

        # this is the only reason it outputs png...
        canvas = FigureCanvasAgg(figure)

        canvas.print_figure(filename)


def plot(f, xmin, xmax, plot_points=25, plot_division=10000000, max_bend=0.1):
    """

    EXAMPLES:
        sage: h1 = lambda x : sqrt(x^3  - 1)
        sage: h2 = lambda x : -sqrt(x^3  - 1)
        sage: plot([h1, h2], 1,5)
    """
    if isinstance(f, (list, tuple)):
        i = 0
        G = plot(f[0], xmin, xmax, plot_points, plot_division, max_bend)
        for i in range(1, len(f)):
            G += plot(f[i], xmin, xmax, plot_points, plot_division, max_bend)
        return G

    xmin = float(xmin)
    xmax = float(xmax)
    plot_points = int(plot_points)
    delta = (xmax - xmin) / plot_points
    xdata = [xmin + i*delta for i in xrange(plot_points)]
    ydata = [float(f(x)) for x in xdata]

    # adaptive refinement
    i = 0
    j = 0
    while i < len(xdata) - 1:
        if abs(ydata[i+1] - ydata[i]) > max_bend:
            x = (xdata[i+1] + xdata[i])/2
            y = float(f(x))
            xdata.insert(i+1, x)
            ydata.insert(i+1, y)
            j += 1
            if j > plot_division:
                break
        else:
            i += 1

    g = Graphics()
    g.line2d(xdata, ydata)

    return g




