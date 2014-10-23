"""
Histograms
"""

#*****************************************************************************
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
from sage.plot.plot import minmax_data, Graphics
from sage.plot.misc import options, rename_keyword

class Histogram(GraphicPrimitive):
    """
    Graphics primitive that represents a histogram.  This takes
    quite a few options as well.

    EXAMPLES::

        sage: from sage.plot.histogram import Histogram
        sage: g = Histogram([1,3,2,0], {}); g
        Histogram defined by a data list of size 4
        sage: type(g)
        <class 'sage.plot.histogram.Histogram'>
        sage: opts = { 'bins':20, 'label':'mydata'}
        sage: g = Histogram([random() for _ in xrange(500)], opts); g
        Histogram defined by a data list of size 500

    We can accept multiple sets of the same length::

        sage: g = Histogram([[1,3,2,0], [4,4,3,3]], {}); g
        Histogram defined by 2 data lists
    """
    def __init__(self, datalist, options):
        """
        Initialize a ``Histogram`` primitive along with
        its options.

        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: Histogram([10,3,5], {'width':0.7})
            Histogram defined by a data list of size 3
        """
        import numpy as np
        self.datalist=np.asarray(datalist,dtype=float)
        GraphicPrimitive.__init__(self, options)        

    def get_minmax_data(self):
        """
        Get minimum and maximum horizontal and vertical ranges
        for the Histogram object.

        EXAMPLES::

            sage: H = histogram([10,3,5], normed=True); h = H[0]
            sage: h.get_minmax_data()
            {'xmax': 10.0, 'xmin': 3.0, 'ymax': 0.4761904761904765, 'ymin': 0}
            sage: G = histogram([random() for _ in xrange(500)]); g = G[0]
            sage: g.get_minmax_data() # random output
            {'xmax': 0.99729312925213209, 'xmin': 0.00013024562219410285, 'ymax': 61, 'ymin': 0}
            sage: Y = histogram([random()*10 for _ in xrange(500)], range=[2,8]); y = Y[0]
            sage: ymm = y.get_minmax_data(); ymm['xmax'], ymm['xmin']
            (8.0, 2.0)
            sage: Z = histogram([[1,3,2,0], [4,4,3,3]]); z = Z[0]
            sage: z.get_minmax_data()
            {'xmax': 4.0, 'xmin': 0, 'ymax': 2, 'ymin': 0}
        """
        import numpy
        options=self.options()
        opt=dict(range = options.pop('range',None),
                 bins = options.pop('bins',None),
                 normed = options.pop('normed',None),
                 weights = options.pop('weights', None))
 
        #check to see if a list of datasets
        if not hasattr(self.datalist[0],'__contains__' ):
            ydata,xdata=numpy.histogram(self.datalist, **opt)
            return minmax_data(xdata,[0]+list(ydata), dict=True)
        else:
            m = { 'xmax': 0, 'xmin':0, 'ymax':0, 'ymin':0}
            for d in self.datalist:
                ydata, xdata = numpy.histogram(d,**opt)
                m['xmax'] = max([m['xmax']] + list(xdata))
                m['xmin'] = min([m['xmin']] + list(xdata))
                m['ymax'] = max([m['ymax']] + list(ydata))
                m['ymin'] = min([m['ymin']] + list(ydata))
            return m
     
    def _allowed_options(self):
        """
        Return the allowed options with descriptions for this graphics
        primitive. This is used in displaying an error message when the
        user gives an option that doesn't make sense.

        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: g = Histogram( [1,3,2,0], {})
            sage: L = list(sorted(g._allowed_options().iteritems()))
            sage: L[0]
            ('align',
             'How the bars align inside of each bin. Acceptable values are "left", "right" or "mid".')
            sage: L[-1]
            ('zorder', 'The layer level to draw the histogram')
            sage: len(L)
            11
        """
        return {'color': 'The color of the face of the bars or list of colors if multiple data sets are given.',
                'edgecolor':'The color of the the border of each bar.',
                'hue':'The color of the bars given as a hue.',
                'zorder':'The layer level to draw the histogram',
                'bins': 'The number of sections in which to divide the range. Also can be a sequence of points within the range that create the partition.', 
                'align': 'How the bars align inside of each bin. Acceptable values are "left", "right" or "mid".',
                'rwidth': 'The relative width of the bars as a fraction of the bin width',
                'range': 'A list [min, max] which define the range of the histogram. Values outside of this range are treated as outliers and omitted from counts.', 
                'normed': '(True or False) If True, the counts are normalized to form a probability density. (n/(len(x)*dbin)', 
                'weights': 'A sequence of weights the same length as the data list. If supplied, then each value contributes it\'s associated weight to the bin count.',
                'label': 'A string label for each data list given.'}

    def _repr_(self):
        """
        Return text representation of this histogram graphics primitive.

        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: g = Histogram( [1,3,2,0], {})
            sage: g._repr_()
            'Histogram defined by a data list of size 4'
            sage: g = Histogram( [[1,1,2,3], [1,3,2,0]], {})
            sage: g._repr_()
            'Histogram defined by 2 data lists'
        """
        L = len(self.datalist)
        if not hasattr(self.datalist[0],'__contains__' ):
            return "Histogram defined by a data list of size {}".format(L)
        else:
            return "Histogram defined by {} data lists".format(L)

    def _render_on_subplot(self, subplot):
        """
        Render this bar chart graphics primitive on a matplotlib subplot
        object.

        EXAMPLES:

        This rendering happens implicitly when the following command
        is executed::

            sage: histogram([1,2,10])  # indirect doctest
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        #check to see if a list of datasets
        if not hasattr(self.datalist[0],'__contains__' ):
            subplot.hist(self.datalist, **options)
        else:
            subplot.hist(self.datalist.transpose(), **options)


@options(aspect_ratio='automatic',align='mid', weights=None, range=None, bins=10, normed=False, edgecolor='black')
def histogram(datalist, **options):
    """
    Computes and draws the histogram for list(s) of numerical data.

    EXAMPLES:

    A very basic histogram for four data points::

        sage: histogram([1,2,3,4], bins=2)
        Graphics object consisting of 1 graphics primitive

    We can see how the histogram compares to various distributions::

        sage: nv = normalvariate
        sage: H = histogram([nv(0,1) for _ in xrange(1000)], bins=20, normed=True, range=[-5,5])
        sage: P = plot( 1/sqrt(2*pi)*e^(-x^2/2), (x,-5,5), color='red', linestyle='--')
        sage: H+P
        Graphics object consisting of 2 graphics primitives

    There are many options one can use with histograms::

        sage: histogram(range(100), color=(1,0,0))
        Graphics object consisting of 1 graphics primitive

    We can do several data sets at once if desired::

        sage: histogram([srange(0,1,.1)*10, [nv(0, 1) for _ in xrange(100)]], color=['red','green'], bins=5)
        Graphics object consisting of 1 graphics primitive
    """

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Histogram(datalist, options=options))
    return g
