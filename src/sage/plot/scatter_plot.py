"""
Scatter Plots
"""

#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Emily Kirkman
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.plot.primitive import GraphicPrimitive
from sage.plot.misc import options

class ScatterPlot(GraphicPrimitive):
    """
    Scatter plot graphics primitive.
    """
    def __init__(self, xdata, ydata, options):
        """
        Scatter plot graphics primitive.

        EXAMPLES::

            sage: import numpy
            sage: from sage.plot.scatter_plot import ScatterPlot
            sage: ScatterPlot(numpy.array([0,1,2]), numpy.array([3.5,2,5.1]), {'facecolor':'white', 'marker':'s'})
            Scatter plot graphics primitive on 3 data points
        """
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: s = scatter_plot([[0,1],[2,4],[3.2,6]])
            sage: d = s.get_minmax_data()
            sage: d['xmin']
            0.0
            sage: d['ymin']
            1.0
        """
        return {'xmin': self.xdata.min(),
                'xmax': self.xdata.max(),
                'ymin': self.ydata.min(),
                'ymax': self.ydata.max()}

    def _allowed_options(self):
        """
        Return the dictionary of allowed options for the scatter plot
        graphics primitive.

        EXAMPLES::

            sage: from sage.plot.scatter_plot import ScatterPlot
            sage: list(sorted(ScatterPlot([-1,2], [17,4], {})._allowed_options().iteritems()))
            [('alpha', 'How transparent the marker border is.'),
             ('edgecolor', 'The color of the marker border.'),
             ('facecolor', 'The color of the marker face.'),
             ('hue', 'The color given as a hue.'),
             ('marker', 'What shape to plot the points.'),
             ('markersize', 'the size of the markers.'),
             ('rgbcolor', 'The color as an rgb tuple.'),
             ('zorder', 'The layer level in which to draw.')]
        """
        return {'markersize': 'the size of the markers.',
                'marker': 'What shape to plot the points.',
                'alpha':'How transparent the marker border is.',
                'rgbcolor':'The color as an rgb tuple.',
                'hue':'The color given as a hue.',
                'facecolor':'The color of the marker face.',
                'edgecolor':'The color of the marker border.',
                'zorder':'The layer level in which to draw.'}

    def _repr_(self):
        """
        Text representation of a scatter plot graphics primitive.

        EXAMPLES::

            sage: import numpy
            sage: from sage.plot.scatter_plot import ScatterPlot
            sage: ScatterPlot(numpy.array([0,1,2]), numpy.array([3.5,2,5.1]), {})
            Scatter plot graphics primitive on 3 data points
        """
        return 'Scatter plot graphics primitive on %s data points'%len(self.xdata)

    def _render_on_subplot(self, subplot):
        """
        Render this scatter plot in a subplot.  This is the key function that
        defines how this scatter plot graphics primitive is rendered in
        matplotlib's library.
        """
        from matplotlib.pyplot import scatter
        options = self.options()
        subplot.scatter(self.xdata, self.ydata, alpha=options['alpha'], zorder=options['zorder'], marker=options['marker'],s=options['markersize'],facecolors=options['facecolor'], edgecolors=options['edgecolor'])

@options(alpha=1, markersize=50, marker='o', zorder=5, facecolor='#fec7b8', edgecolor='black')
def scatter_plot(datalist, **options):
    """
    Returns a Graphics object of a scatter plot containing all points in
    the datalist.  Type \code{scatter_plot.options} to see all available
    plotting options.

    INPUT:

    - ``datalist`` -- a list of tuples ``(x,y)``

    - ``alpha`` -- default: 1

    - ``markersize`` -- default: 50

    - ``marker`` -- default: ``'o'``

    - ``facecolor`` -- default: ``'#fec7b8'``

    - ``edgecolor`` -- default: ``'black'``

    - ``zorder`` -- default: 5

    EXAMPLES::

        sage: s = scatter_plot([[0,1],[2,2],[4.3,1.1]], marker='s')
        sage: s

    """
    import numpy
    from sage.plot.plot import Graphics
    g = Graphics()
    data = numpy.array(datalist, dtype='float')
    if len(data) != 0:
        xdata = data[:,0]
        ydata = data[:,1]
        g.add_primitive(ScatterPlot(xdata, ydata, options=options))
    return g
