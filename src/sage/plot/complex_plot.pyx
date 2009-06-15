"""
Complex Plots
"""

#*****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>,
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

# TODO: use NumPy buffers and complex fast_callable (when supported)

include "../ext/stdsage.pxi"
include "../ext/interrupt.pxi"

cimport numpy as cnumpy
import numpy

from sage.plot.primitive import GraphicPrimitive
from sage.plot.misc import options
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.misc.misc import srange

cdef extern from "math.h":
    double hypot(double, double)
    double atan2(double, double)
    double atan(double)
    double log(double)
    double exp(double)
    double PI

cdef inline ComplexDoubleElement new_CDF_element(double x, double y):
    cdef ComplexDoubleElement z = <ComplexDoubleElement>PY_NEW(ComplexDoubleElement)
    z._complex.dat[0] = x
    z._complex.dat[1] = y
    return z

cdef inline double mag_to_lightness(double r):
    """
    Tweak this to adjust how the magnitude affects the color.

    INPUT:

    - ``r`` - a non-negative real number

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.
    """
    return atan(log(r+1)) * (4/PI) - 1

def complex_to_rgb(z_values):
    """
    INPUT:

    - ``z_values`` -- A grid of complex numbers, as a list of lists

    OUTPUT:

    An `N \\times M \\times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r,g,b) tuple.

    EXAMPLES::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 1000]])
        array([[[ 0.        ,  0.        ,  0.        ],
                [ 0.77172568,  0.        ,  0.        ],
                [ 1.        ,  0.81697746,  0.81697746]]])
        sage: complex_to_rgb([[0, 1j, 1000j]])
        array([[[ 0.        ,  0.        ,  0.        ],
                [ 0.38586284,  0.        ,  0.77172568],
                [ 0.90848873,  0.81697746,  1.        ]]])
    """
    import numpy
    cdef unsigned int i, j, imax, jmax
    cdef double x, y, mag, arg
    cdef double lightness, hue, top, bot
    cdef double r, g, b
    cdef int ihue
    cdef ComplexDoubleElement z
    from sage.rings.complex_double import CDF

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] rgb = numpy.empty(dtype=numpy.float, shape=(imax, jmax, 3))

    _sig_on
    for i from 0 <= i < imax:

        row = z_values[i]
        for j from 0 <= j < jmax:

            zz = row[j]
            z = <ComplexDoubleElement>(zz if PY_TYPE_CHECK_EXACT(zz, ComplexDoubleElement) else CDF(zz))
            x, y = z._complex.dat[0], z._complex.dat[1]
            mag = hypot(x, y)
            arg = atan2(y, x)

            lightness = mag_to_lightness(mag)
            if lightness < 0:
                bot = 0
                top = (1+lightness)
            else:
                bot = lightness
                top = 1

            hue = -3*arg/PI
            if hue < 0: hue += 6
            ihue = <int>hue
            if ihue == 0:
                r = top
                g = bot + hue * (top-bot)
                b = bot
            elif ihue == 1:
                r = bot + (2-hue) * (top-bot)
                g = top
                b = bot
            elif ihue == 2:
                r = bot
                g = top
                b = bot + (hue-2) * (top-bot)
            elif ihue == 3:
                r = bot
                g = bot + (4-hue) * (top-bot)
                b = top
            elif ihue == 4:
                r = bot + (hue-4) * (top-bot)
                g = bot
                b = top
            else:
                r = top
                g = bot
                b = bot + (6-hue) * (top-bot)

            rgb[i, j, 0] = r
            rgb[i, j, 1] = g
            rgb[i, j, 2] = b

    _sig_off
    return rgb


class ComplexPlot(GraphicPrimitive):
    def __init__(self, z_values, xrange, yrange, options):
        """
        TESTS::

            sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
        """
        self.xrange = xrange
        self.yrange = yrange
        self.z_values = z_values
        self.x_count = len(z_values)
        self.y_count = len(z_values[0])
        self.rgb_data = complex_to_rgb(z_values)
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: p = complex_plot(lambda z: z, (-1, 2), (-3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', -1.0), ('ymax', 4.0), ('ymin', -3.0)]
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xrange, self.yrange, dict=True)

    def _allowed_options(self):
        """
        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._allowed_options(), dict)
            True
        """
        return {'plot_points':'How many points to use for plotting precision',
                'interpolation':'What interpolation method to use'}

    def _repr_(self):
        """
        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._repr_(), str)
            True
        """
        return "ComplexPlot defined by a %s x %s data grid"%(self.x_count, self.y_count)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: complex_plot(lambda x: x^2, (-5, 5), (-5, 5))
        """
        options = self.options()
        x0,x1 = float(self.xrange[0]), float(self.xrange[1])
        y0,y1 = float(self.yrange[0]), float(self.yrange[1])
        subplot.imshow(self.rgb_data, origin='lower', extent=(x0,x1,y0,y1), interpolation=options['interpolation'])

@options(plot_points=100, interpolation='catrom')
def complex_plot(f, xrange, yrange, **options):
    r"""
    ``complex_plot`` takes a complex function of one variable,
    `f(z)` and plots output of the function over the specified
    ``xrange`` and ``yrange`` as demonstrated below. The magnitude of the
    output is indicated by the brightness (with zero being black and
    infinity being white) while the argument is represented by the
    hue (with red being positive real, and increasing through orange,
    yellow, ... as the argument increases.

    ``complex_plot(f, (xmin, xmax), (ymin, ymax), ...)``

    INPUT:

    - ``f`` -- a function of a single complex value `x + iy`

    - ``(xmin, xmax)`` -- 2-tuple, the range of ``x`` values

    - ``(ymin, ymax)`` -- 2-tuple, the range of ``y`` values

    The following inputs must all be passed in as named parameters:

    - ``plot_points`` -- integer (default: 100); number of points to plot
      in each direction of the grid

    - ``interpolation`` -- string (default: ``'catrom'``), the interpolation
      method to use: ``'bilinear'``, ``'bicubic'``, ``'spline16'``,
      ``'spline36'``, ``'quadric'``, ``'gaussian'``, ``'sinc'``,
      ``'bessel'``, ``'mitchell'``, ``'lanczos'``, ``'catrom'``,
      ``'hermite'``, ``'hanning'``, ``'hamming'``, ``'kaiser'``


    EXAMPLES:

    Here we plot a couple of simple functions::

        sage: complex_plot(sqrt, (-5, 5), (-5, 5))
        sage: complex_plot(sin, (-5, 5), (-5, 5))
        sage: complex_plot(log, (-10, 10), (-10, 10))
        sage: complex_plot(exp, (-10, 10), (-10, 10))

    A function with some nice zeros and a pole::

        sage: f(z) = z^5 + z - 1 + 1/z
        sage: complex_plot(f, (-3, 3), (-3, 3))

    Here is the identity, useful for seeing what values map to what colors::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3))

    The Riemann Zeta function::

        sage: complex_plot(zeta, (-30,30), (-30,30))
    """
    from sage.plot.plot import Graphics, setup_for_eval_on_grid
    cdef double x, y
    ignore, xstep, ystep, xrange, yrange = setup_for_eval_on_grid([f], xrange, yrange, options['plot_points'])
    xmin, xmax = xrange
    ymin, ymax = yrange
    xrange_list = srange(xmin, xmax+xstep, xstep, universe=float)
    yrange_list = srange(ymin, ymax+ystep, ystep, universe=float)
    _sig_on
    z_values = [[  f(new_CDF_element(x, y)) for x in xrange_list]
                                            for y in yrange_list]
    _sig_off
    g = Graphics()
    g.add_primitive(ComplexPlot(z_values, xrange, yrange, options))
    return g
