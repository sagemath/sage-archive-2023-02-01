"""
Histogram Charts
"""

#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>, 
#                     2010 Jason Grout <jason-sage@creativetrax.com>
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
from sage.plot.plot import minmax_data, Graphics
from sage.plot.misc import options, rename_keyword

class Histogram(GraphicPrimitive):
    """
    Graphics primitive that represents a histogram.

    EXAMPLES::

        sage: from sage.plot.bar_chart import BarChart
        sage: g = Histogram(range(4), [1,3,2,0], {}); g
        Histogram defined by a 4 datalist 
        sage: type(g)
        <class 'sage.plot.histogram.Histogram'>
    """
    def __init__(self, datalist, options):
        """
        Initialize a ``Histogram`` primitive.
        
        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: Histogram(range(3), [10,3,5], {'width':0.7})
            Histogram defined by a 3 datalist 
        """
        self.datalist = datalist
        GraphicPrimitive.__init__(self, options)        

    def get_minmax_data(self):
        """
        """
        import numpy
        options=self.options()
        opt=dict(range=options['range'],
                 bins=options['bins'],
                 normed=options['normed'],
                 weights=options['weights'])
 
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
            sage: g = Histogram(range(4), [1,3,2,0], {})
            sage: list(sorted(g._allowed_options().iteritems()))
            [('hue', 'The color given as a hue.'),
             ('rgbcolor', 'The color as an RGB tuple.'),
             ('width', 'The width of the bars'),
             ('zorder', 'The layer level in which to draw')]
        """
        return {'color':'The color as an RGB tuple.',
                'hue':'The color given as a hue.',
                'zorder':'The layer level in which to draw',
                'bins': 'asdf', 'align': 'asf', 'range': 'asf', 'normed': 'asf', 'weights': 'asf','edgecolor':'asdf'}

    def _repr_(self):
        """
        Return text representation of this histogram graphics primitive.

        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: g = Histogram(range(4), [1,3,2,0], {})
            sage: g._repr_()
            'BarChart defined by a 4 datalist'         
        """
        return "BarChart defined by a %s datalist"%(len(self.datalist))

    def _render_on_subplot(self, subplot):
        """
        Render this bar chart graphics primitive on a matplotlib subplot
        object.

        EXAMPLES:

        This rendering happens implicitly when the following command
        is executed::

            sage: bar_chart([1,2,10])
        """
        options = self.options()
        subplot.hist(self.datalist, **options)

@rename_keyword(rgbcolor='color')
@options(color='lightblue',align='mid', weights=None, range=None, normed=False, bins=50,edgecolor='black')
def histogram(datalist, **options):
    """
    A bar chart of (currently) one list of numerical data.
    Support for more data lists in progress.

    EXAMPLES:

    A bar_chart with blue bars::

        sage: bar_chart([1,2,3,4])

    A bar_chart with thinner bars::

        sage: bar_chart([x^2 for x in range(1,20)], width=0.2)

    A bar_chart with negative values and red bars::

        sage: bar_chart([-3,5,-6,11], rgbcolor=(1,0,0))

    Extra options will get passed on to show(), as long as they are valid::

        sage: bar_chart([-2,8,-7,3], rgbcolor=(1,0,0), axes=False)
        sage: bar_chart([-2,8,-7,3], rgbcolor=(1,0,0)).show(axes=False) # These are equivalent
    """

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Histogram(datalist, options=options))
    return g
