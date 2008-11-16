#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
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
from sage.plot.misc import options, rename_keyword, to_mpl_color
from sage.misc.misc import xsrange

# Below is the base class that is used to make 'field plots'.
# Its implementation is motivated by 'PlotField'.
# Currently it is used to make the function 'plot_vector_field'
# TODO: use this to make these functions:
# 'plot_gradient_field' and 'plot_hamiltonian_field'
class PlotField(GraphicPrimitive):
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

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES:
            sage: x,y = var('x,y')
            sage: d = plot_vector_field((.01*x,x+y), (10,20), (10,20))[0].get_minmax_data()
            sage: d['xmin']
            10.0
            sage: d['ymin']
            10.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xpos_array, self.ypos_array, dict=True)

    def _allowed_options(self):
        return {'plot_points':'How many points to use for plotting precision',
                'pivot': 'Where the arrow should be placed in relation to the point (tail, middle, tip)',
                'headwidth': 'Head width as multiple of shaft width, default is 3',
                'headlength': 'head length as multiple of shaft width, default is 5',
                'headaxislength': 'head length at shaft intersection, default is 4.5'}

    def _repr_(self):
        return "PlotField defined by a %s x %s vector grid"%(len(self.xpos_array), len(self.ypos_array))

    def _render_on_subplot(self, subplot):
        options = self.options()
        quiver_options = options.copy()
        quiver_options.pop('plot_points')
        subplot.quiver(self.xpos_array, self.ypos_array, self.xvec_array, self.yvec_array, angles='xy', **quiver_options)

@options(plot_points=20)
def plot_vector_field((f, g), xrange, yrange, **options):
    r"""

    \code{plot_vector_field} takes two functions of two variables, $(f(x,y), g(x,y))$
    and plots vector arrows of the function over the specified
    xrange and yrange as demonstrated below.

    plot_vector_field((f, g), (xvar, xmin, xmax), (yvar, ymin, ymax))

    EXAMPLES:
    Plot the vector fields involving sin and cos
        sage: x,y = var('x y')
        sage: plot_vector_field((sin(x), cos(y)), (x,-3,3), (y,-3,3))
        sage: plot_vector_field(( y, (cos(x)-2)*sin(x)), (x,-pi,pi), (y,-pi,pi))

    Plot a gradient field
        sage: u,v = var('u v')
        sage: f = exp(-(u^2+v^2))
        sage: plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2))


    """
    from sage.plot.plot import setup_for_eval_on_grid, Graphics
    z, xstep, ystep, xrange, yrange = setup_for_eval_on_grid([f,g], xrange, yrange, options['plot_points'])
    f,g = z

    xpos_array, ypos_array, xvec_array, yvec_array = [],[],[],[]
    for x in xsrange(xrange[0], xrange[1], xstep):
        for y in xsrange(yrange[0], yrange[1], ystep):
            xpos_array.append(x)
            ypos_array.append(y)
            xvec_array.append(f(x,y))
            yvec_array.append(g(x,y))

    import numpy
    xvec_array = numpy.array(xvec_array, dtype=float)
    yvec_array = numpy.array(yvec_array, dtype=float)
    g = Graphics()
    g.add_primitive(PlotField(xpos_array, ypos_array, xvec_array, yvec_array, options))
    return g

def plot_slope_field(f, xrange, yrange, **kwds):
    r"""

    \code{plot_slope_field} takes a function of two variables, $f(x,y)$, and at various points (x_i,y_i), plots a line with slope $f(x_i,y_i)$

    plot_slope_field((f, g), (xvar, xmin, xmax), (yvar, ymin, ymax))

    EXAMPLES:
    A logistic function modeling population growth.
        sage: x,y = var('x y')
        sage: capacity = 3 # thousand
        sage: growth_rate = 0.7 # population increases by 70% per unit of time
        sage: plot_slope_field(growth_rate*(1-y/capacity)*y, (x,0,5), (y,0,capacity*2)).show(aspect_ratio=1)

    Plot a slope field involving sin and cos
        sage: x,y = var('x y')
        sage: plot_slope_field(sin(x+y)+cos(x+y), (x,-3,3), (y,-3,3)).show(aspect_ratio=1)

    """
    from math import sqrt
    slope_options = {'headaxislength': 0, 'headlength': 0, 'pivot': 'middle'}
    slope_options.update(kwds)

    from sage.calculus.calculus import sqrt
    norm = sqrt((f**2+1))
    return plot_vector_field((1/norm, f/norm), xrange, yrange, **slope_options)
