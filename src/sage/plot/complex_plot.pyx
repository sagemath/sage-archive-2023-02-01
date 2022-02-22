"""
Complex Plots

AUTHORS:

- Robert Bradshaw (2009): initial version
- David Lowry-Duda (2022): incorporate matplotlib colormaps
"""
# ****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>,
#       Copyright (C) 2022 David Lowry-Duda <david@lowryduda.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# TODO: use NumPy buffers and complex fast_callable (when supported)
from cysignals.signals cimport sig_on, sig_off, sig_check

# XXX: Temporary import, remove when cythonized
import colorsys

cimport numpy as cnumpy

from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.arith.srange import srange

from libc.math cimport hypot, atan2, atan, log, log2, sqrt
from sage.arith.constants cimport M_PI as PI

from sage.libs.gsl.complex cimport *


cdef inline ComplexDoubleElement new_CDF_element(double x, double y):
    z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
    GSL_SET_COMPLEX(&z._complex, x, y)
    return z


cdef inline double mag_to_lightness(double r):
    """
    Returns a lightness for the given magnitude.

    Small magnitudes are darker and large magnitudes are lighter, and the
    lightness smoothly transitions.

    Tweak this to adjust how the magnitude affects the color. For instance,
    changing ``sqrt(r)`` to ``r`` will cause anything near a zero to be much
    darker and poles to be much lighter, while ``r**(.25)`` would cause the
    reverse effect.

    INPUT:

    - ``r`` - a non-negative real number

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.cyclic_mag_to_lightness`
        :func:`sage.plot.complex_plot.mag_and_arg_to_lightness`

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import direct_complex_to_rgb
        sage: direct_complex_to_rgb([[0, 1, 10]])
        array([[[0.        , 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ],
                [1.        , 0.22134776, 0.22134776]]])
    """
    return atan(log(sqrt(r)+1)) * (4/PI) - 1


cdef inline double cyclic_mag_to_lightness(double r):
    """
    Returns a lightness for the given magnitude.

    This modifies the lightness around magnitudes of size `2^n` for integer `n`
    to create the appearance of logarithmically spaced contours.

    Tweaking this will affect how the magnitude affects the color. Changing
    ``log2(r)`` to ``log2(r)/log2(5)`` will associate contours around
    magnitudes of size `5^n` for integer `n`. Changing ``log2(r)`` to ``r``
    will give cyclic linear contours, which might be more appropriate for
    complex functions with very little variation in size.

    INPUT:

    - ``r`` - a non-negative real number

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.mag_to_lightness`
        :func:`sage.plot.complex_plot.mag_and_arg_to_lightness`

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import direct_complex_to_rgb
        sage: direct_complex_to_rgb([[0, 1, 10]], contoured=True)
        array([[[0.65      , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [0.98903595, 0.        , 0.        ]]])
    """
    if r < 1e-10:
        return -0.35
    rem = log2(r) % 1
    if rem < 0:  # Choose positive mod representative
        rem += 1
    return .15 - rem/2.


cdef inline double mag_and_arg_to_lightness(double r, double arg):
    """
    Returns a lightness for the given magnitude and argument.

    This modifies the lightness around magnitudes of size `2^n` for integer
    `n`, and also around arguments of the form `(2 \pi) / 10 * n` for integer
    `n`. This creates the appearance of tiles, bounded by logarithmically
    spaced magnitude contours and evenly spaced argument contours.

    Tweaking the magnitude function ``log2(r)` will affect the spacing of
    magnitude contours. Tweaking the argument computation, such as ``(5 * arg /
    PI)`` to ``(2 * arg / PI)`` will change the number of phase contours from
    10 to 4.

    INPUT:

    - ``r`` - a non-negative real number
    - ``arg`` - a real number

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.mag_to_lightness`
        :func:`sage.plot.complex_plot.cyclic_mag_to_lightness`

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import direct_complex_to_rgb
        sage: direct_complex_to_rgb([[0, 1, 10]], tiled=True)
        array([[[0.65      , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [1.        , 0.06951798, 0.06951798]]])
    """
    if r < 1e-10:
        return -0.35
    cdef double r_rem, arg_rem
    r_rem = log2(r) % 1
    arg_rem = (5 * arg / PI) % 1
    if r_rem < 0:  # Choose positive mod representatives
        r_rem += 1
    if arg_rem < 0:
        arg_rem += 1
    return 0.15 - (r_rem)/4. - (arg_rem)/4.


def direct_complex_to_rgb(z_values, contoured=False, tiled=False):
    """
    INPUT:

    - ``z_values`` -- A grid of complex numbers, as a list of lists

    - ``contoured`` -- boolean (default: False) - causes magnitude to be
      indicated through contour-like adjustments to lightness.

    - ``tiled`` -- boolean (default: False) - causes magnitude and argument to
      be indicated through contour-like adjustments to lightness.

    OUTPUT:

    An `N \\times M \\times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r,g,b) tuple.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.cmap_complex_to_rgb`

    EXAMPLES::

        sage: from sage.plot.complex_plot import direct_complex_to_rgb
        sage: direct_complex_to_rgb([[0, 1, 1000]])
        array([[[0.        , 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ],
                [1.        , 0.64421177, 0.64421177]]])
        sage: direct_complex_to_rgb([[0, 1j, 1000j]])
        array([[[0.        , 0.        , 0.        ],
                [0.38586284, 0.77172568, 0.        ],
                [0.82210588, 1.        , 0.64421177]]])
        sage: direct_complex_to_rgb([[0, 1, 1000]], contoured=True)
        array([[[0.65      , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [0.66710786, 0.        , 0.        ]]])
        sage: direct_complex_to_rgb([[0, 1, 1000]], tiled=True)
        array([[[0.65      , 0.        , 0.        ],
                [1.        , 0.15      , 0.15      ],
                [0.90855393, 0.        , 0.        ]]])
    """
    import numpy as np

    cdef unsigned int i, j, imax, jmax
    cdef double x, y, mag, arg
    cdef double lightness, hue, top, bot
    cdef double r, g, b
    cdef int ihue
    cdef ComplexDoubleElement z
    from sage.rings.complex_double import CDF

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] rgb = np.empty(dtype=float, shape=(imax, jmax, 3))

    sig_on()
    for i from 0 <= i < imax:

        row = z_values[i]
        for j from 0 <= j < jmax:

            zz = row[j]
            if type(zz) is ComplexDoubleElement:
                z = <ComplexDoubleElement>zz
            else:
                z = CDF(zz)
            x = GSL_REAL(z._complex)
            y = GSL_IMAG(z._complex)
            mag = hypot(x, y)
            arg = atan2(y, x) # math module arctan has range from -pi to pi, so cut along negative x-axis

            if tiled:
                lightness = mag_and_arg_to_lightness(mag, arg)
            elif contoured:
                lightness = cyclic_mag_to_lightness(mag)
            else:
                lightness = mag_to_lightness(mag)

            if lightness < 0: # in hsv, variable value, full saturation (s=1, v=1+lightness)
                bot = 0
                top = (1+lightness)
            else: # in hsv, variable saturation, full value (v=1, s=1-lightness)
                bot = lightness
                top = 1

            # Note that does same thing as colorsys module hsv_to_rgb
            # for this setup, but in Cython
            hue = 3*arg/PI

            if hue < 0: hue += 6 # usual hsv hue is thus h=arg/(2*pi) for positive, h=arg/(2*PI)+1 for negative
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

    sig_off()
    return rgb


def cmap_complex_to_rgb(z_values, cmap=None, contoured=False, tiled=False):
    """
    INPUT:

    - ``z_values`` -- A grid of complex numbers, as a list of lists

    - ``contoured`` -- boolean (default: False) - causes magnitude to be
      indicated through contour-like adjustments to lightness.

    - ``tiled`` -- boolean (default: False) - causes magnitude and argument to
      be indicated through contour-like adjustments to lightness.

    OUTPUT:

    An `N \\times M \\times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r, g, b) tuple.

    .. SEEALSO::

        :func:`sage.plot.complex_plot.direct_complex_to_rgb`

    EXAMPLES:

        sage: import matplotlib.cm
        sage: from sage.plot.complex_plot import cmap_complex_to_rgb
        sage: cmap_complex_to_rgb([[0, 1, 1000]])
        array([[[0.        , 0.        , 0.        ],
                [0.49669808, 0.76400071, 0.18024425],
                [0.87320419, 0.99643856, 0.72730967]]])
        sage: cmap_complex_to_rgb([[0, 1, 1000]], cmap=matplotlib.cm.viridis)
        array([[[0.        , 0.        , 0.        ],
                [0.0984475 , 0.4375291 , 0.42487821],
                [0.68959896, 0.84592555, 0.84009311]]])
    """
    import numpy as np
    import matplotlib.cm
    if cmap is None or cmap == 'default':
        cmap = matplotlib.cm.get_cmap('turbo')

    cdef unsigned int i, j, imax, jmax
    cdef double x, y, mag, arg
    cdef double lightness_delta
    cdef double r, g, b
    cdef ComplexDoubleElement z
    from sage.rings.complex_double import CDF

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] als = np.empty(dtype=float, shape=(imax, jmax, 2))
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] rgbs = np.empty(dtype=float, shape=(imax, jmax, 3))

    sig_on()
    for i from 0 <= i < imax:
        row = z_values[i]
        for j from 0 <= j < jmax:
            zz = row[j]
            if type(zz) is ComplexDoubleElement:
                z = <ComplexDoubleElement>zz
            else:
                z = CDF(zz)
            x = GSL_REAL(z._complex)
            y = GSL_IMAG(z._complex)
            mag = hypot(x, y)
            arg = atan2(y, x) # math module arctan has range from -pi to pi, so cut along negative x-axis
            if tiled:
                lightness_delta = mag_and_arg_to_lightness(mag, arg)
            elif contoured:
                lightness_delta = cyclic_mag_to_lightness(mag)
            else:
                # Note in this case that `lightness_delta` is interpreted
                # differently by `set_darkness` in comparison to above.
                lightness_delta = mag_to_lightness(mag)
            als[i, j, 0] = arg
            als[i, j, 1] = lightness_delta

    args = als[:,:,0]
    nan_indices = np.isnan(als).any(-1)            # Mask for undefined points
    normalized_colors = cmap((args + PI) / (2 * PI))
    normalized_colors = normalized_colors[:,:,:3]  # discard alpha channel
    lightdeltas = als[:,:,1]
    lightdeltas = lightdeltas # normalize lightness adjustment
    arg_d_s = np.dstack((normalized_colors, lightdeltas))

    if tiled or contoured:
        # add contours: convert to hls, adjust lightness, convert back
        tmparr = np.apply_along_axis(rgbld_to_hlsld, 2, arg_d_s)
        tmparr = np.apply_along_axis(adjust_lightness, 2, tmparr)
        rgbs = np.apply_along_axis(hls_to_rgb, 2, tmparr)
        # XXX: These calls are very slow! Implement the following
        #      directly in cython instead.
        # rgbs = manual_contoured_cmap_complex_to_rgb(arg_d_s)
    else:
        rgbs = np.apply_along_axis(set_darkness, 2, arg_d_s)
        # XXX: These calls are very slow! Implement the following
        #      directly in cython instead.
        # rgbs = manual_smooth_cmap_complex_to_rgb(arg_d_s)

    # Apply mask, making nan_indices white
    rgbs[nan_indices] = 1

    sig_off()
    return rgbs


# XXX: This is a python function, and being repeatedly called by
#      np.apply_along_axis is slow. This is to be replaced.
def adjust_lightness(hlsd_arr):
    """
    Create contours by adjusting lightness.

    This applies the lightdelta to the lightness and returns the modified
    `(h, l, s)` as a Numpy array.

    INPUT:

    - ``hlsd_arr`` - a length 4 array `(h, l, s, d)` containing a hue,
      lightness, saturation, and lightnessdelta. Each of `h`, `l`, and `s` are
      between `0` and `1` (inclusive), and `d` is between `-1` and `1`
      (inclusive).

    OUTPUT:

    A Numpy array `(h, l, s)` with each coordinate between `0` and `1`
    (inclusive).

    EXAMPLES:

        sage: from sage.plot.complex_plot import adjust_lightness
        sage: adjust_lightness([0, 0.25, 0.5, 0.75])
        array([0.   , 0.625, 0.5  ])
        sage: adjust_lightness([0, 0.25, 0.5, 1])
        array([0.  , 0.75, 0.5 ])
    """
    import numpy as np
    h, l, s, delta = hlsd_arr
    lightdelta = 0.5*(delta)
    l += lightdelta
    if l < 0:
        l = 0
    if l > 1:
        l = 1
    return np.array([h, l, s])


# XXX: This is a python function, and being repeatedly called by
#      np.apply_along_axis is slow. This is to be replaced.
def set_darkness(rgbmag_arr):
    """
    Adjusts an rgb color darker or lighter based on given magnitude.

    This produces similar gradation to the default (no colormap, no contours)
    complex_plot, but is intended to apply to plots with colormaps.

    INPUT:

    - ``rgbmag_arr`` - a length 4 array `(r, g, b, mag)` containing the rgb
      data and associated magnitude. Each of `r`, `g`, and `g` are between `0`
      and `1` (inclusive). We assume `mag` has been mapped to the range `-1` to
      `1`.

    OUTPUT:

    A Numpy array `(r, g, b)` with each coordinate between `0` and `1`
    (inclusive).

    .. SEEALSO::

        :func:`sage.plot.complex_plot.mag_to_lightness`

    ALGORITHM:

    This is naive. First the magnitude is remapped to the range `[0, 1]`. Then
    this piece-wise linearly interpolates between the given `(r, g, b)` and
    black (when the mapped magnitude is `0`) and white (when the mapped
    magnitude is `1`).

    EXAMPLES:

        sage: from sage.plot.complex_plot import set_darkness
        sage: set_darkness([0, 0.25, 0.5, 0.75])
        array([0.75  , 0.8125, 0.875 ])
        sage: set_darkness([0, 0.25, 0.5, 1])
        array([1., 1., 1.])
        sage: set_darkness([0, 0.25, 0.5, 0.25])
        array([0.25  , 0.4375, 0.625 ])
    """
    import numpy as np

    r, g, b, mag = rgbmag_arr
    mag = 0.5 + mag/2.
    if mag < 0.5:
        return 2*mag*np.array([r, g, b])
    white = np.array([1, 1, 1])
    return 2*(mag - 1)*(white - np.array([r, g, b])) + white


# XXX: This is a python function, and being repeatedly called by
#      np.apply_along_axis is slow. This is to be replaced.
def rgbld_to_hlsld(rgbld_arr):
    """
    Converts `(r, g, b, ld)` to `(h, l, s, ld)`.

    INPUT:

    - ``rgbld_arr`` -- An array `(r, g, b, ld)` of rgb colordata and a
      lightness delta. Each of `r`, `g`, and `b` are between `0` and `1`
      (inclusive), and `ld` is between `-1` and `1` (inclusive).

    OUTPUT:

    A Numpy array `[h, l, s, ld]`, where each of `h`, `l`, and `s` are between
    `0` and `1` (inclusive), and `ld` is between `-1` and `1` (inclusive).

    EXAMPLES:

        sage: from sage.plot.complex_plot import rgbld_to_hlsld
        sage: rgbld_to_hlsld([0, 0.25, 0.5, 0.75])
        array([0.58333333, 0.25      , 1.        , 0.75      ])
    """
    import numpy as np
    r = float(rgbld_arr[0])
    g = float(rgbld_arr[1])
    b = float(rgbld_arr[2])
    ld = float(rgbld_arr[3])
    return np.array(list(colorsys.rgb_to_hls(r, g, b)) + [ld])


# XXX: This is a python function, and being repeatedly called by
#      np.apply_along_axis is slow. This is to be replaced.
def hls_to_rgb(hls_arr):
    """
    Converts `(h, l, s`) to `(r, g, b)`.

    INPUT:

    - ``hls_arr` -- An array `(h, l, s)` of hue, lightness, saturation, each
      between `0` and `1` (inclusive).

    OUTPUT:

    A Numpy array `[r, g, b]`, where each coordinate is between `0` and `1`
    (inclusive).

    EXAMPLES:

        sage: from sage.plot.complex_plot import hls_to_rgb
        sage: hls_to_rgb([0, 0.25, 0.5])
        array([0.375, 0.125, 0.125])

    .. WARNING::

        This behaves differently than ``colorsys.hls_to_rgb`` called directly
        from sage, but **not** differently than ``colorsys.hls_to_rgb`` called
        from a python interpreter. This is due to a difference in rounding
        behavior, for example whether ``0.75 % 1`` is `0.75` or `-0.25`.
    """
    import numpy as np
    cdef double h, l, s
    h = hls_arr[0]
    l = hls_arr[1]
    s = hls_arr[2]
    return np.array(colorsys.hls_to_rgb(h, l, s))


class ComplexPlot(GraphicPrimitive):
    """
    The GraphicsPrimitive to display complex functions in using the domain
    coloring method

    INPUT:

    - ``rgb_data`` -- An array of colored points to be plotted.

    - ``x_range`` -- A minimum and maximum x value for the plot.

    - ``y_range`` -- A minimum and maximum y value for the plot.

    TESTS::

        sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
    """
    def __init__(self, rgb_data, x_range, y_range, options):
        """
        TESTS::

            sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
        """
        self.x_range = x_range
        self.y_range = y_range
        self.x_count = len(rgb_data)
        self.y_count = len(rgb_data[0])
        self.rgb_data = rgb_data
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: p = complex_plot(lambda z: z, (-1, 2), (-3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', -1.0), ('ymax', 4.0), ('ymin', -3.0)]
            sage: p = complex_plot(lambda z: z, (1, 2), (3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', 1.0), ('ymax', 4.0), ('ymin', 3.0)]
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.x_range, self.y_range, dict=True)

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
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        x0, x1 = float(self.x_range[0]), float(self.x_range[1])
        y0, y1 = float(self.y_range[0]), float(self.y_range[1])
        subplot.imshow(self.rgb_data, origin='lower', extent=(x0, x1, y0, y1),
                       interpolation=options['interpolation'])


@options(plot_points=100, interpolation='catrom')
def complex_plot(f, x_range, y_range, contoured=False, tiled=False, cmap=None, **options):
    r"""
    ``complex_plot`` takes a complex function of one variable,
    `f(z)` and plots output of the function over the specified
    ``x_range`` and ``y_range`` as demonstrated below. The magnitude of
    the output is indicated by the brightness and the argument is
    represented by the hue.

    By default, zero magnitude corresponds to black output, infinite
    magnitude corresponds to white output. The options ``contoured``,
    ``tiled``, and ``cmap`` affect the output.

    ``complex_plot(f, (xmin, xmax), (ymin, ymax), contoured, tiled, cmap, ...)``

    INPUT:

    - ``f`` -- a function of a single complex value `x + iy`

    - ``(xmin, xmax)`` -- 2-tuple, the range of ``x`` values

    - ``(ymin, ymax)`` -- 2-tuple, the range of ``y`` values

    - ``contoured`` -- boolean (default: False) - causes the magnitude
      to be indicated by logarithmically spaced 'contours'. The
      magnitude along one contour is either twice or half the magnitude
      along adjacent contours.

    - ``tiled`` -- boolean (default: False) - causes the magnitude to
      be indicated by logarithmically spaced 'contours' as in
      ``contoured``, and in addition for there to be `10` evenly
      spaced phase contours.

    - ``cmap`` --  the string name of a matplotlib colormap, or an instance of
      a matplotlib Colormap, or the special string `default` (default: None).
      If None, then hues are chosen from a standard color wheel, cycling from
      red to yellow to blue. If `default`, then a default matplotlib colormap
      is chosen.

    The following inputs may be passed in as named parameters:

    - ``plot_points`` -- integer (default: 100); number of points to
      plot in each direction of the grid

    - ``interpolation`` -- string (default: ``'catrom'``), the interpolation
      method to use: ``'bilinear'``, ``'bicubic'``, ``'spline16'``,
      ``'spline36'``, ``'quadric'``, ``'gaussian'``, ``'sinc'``,
      ``'bessel'``, ``'mitchell'``, ``'lanczos'``, ``'catrom'``,
      ``'hermite'``, ``'hanning'``, ``'hamming'``, ``'kaiser'``

    .. NOTE::

        Matplotlib colormaps can be chosen or customized to cater to different
        types of vision. The colormaps 'cividis' and 'viridis' in matplotlib
        are designed to be perceptually uniform to a broader audience. The
        colormap 'turbo' is similar to the default but with more even contrast.
        See [NAR2018]_ for more information about colormap choice for
        scientific visualization.


    EXAMPLES:

    Here we plot a couple of simple functions::

        sage: complex_plot(sqrt(x), (-5, 5), (-5, 5))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(sqrt(x), (-5, 5), (-5, 5)))

    ::

        sage: complex_plot(sin(x), (-5, 5), (-5, 5))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(sin(x), (-5, 5), (-5, 5)))

    ::

        sage: complex_plot(log(x), (-10, 10), (-10, 10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(log(x), (-10, 10), (-10, 10)))

    ::

        sage: complex_plot(exp(x), (-10, 10), (-10, 10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(exp(x), (-10, 10), (-10, 10)))

    A plot with a different choice of colormap.

        sage: complex_plot(exp(x), (-10, 10), (-10, 10), cmap='viridis')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(exp(x), (-10, 10), (-10, 10), cmap='viridis))

    A function with some nice zeros and a pole::

        sage: f(z) = z^5 + z - 1 + 1/z
        sage: complex_plot(f, (-3, 3), (-3, 3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3)))

    The same function as above, but with contours.

        sage: f(z) = z^5 + z - 1 + 1/z
        sage: complex_plot(f, (-3, 3), (-3, 3), contoured=True)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3), contours=True))

    Here is the identity, useful for seeing what values map to what colors::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(lambda z: z, (-3, 3), (-3, 3)))

    The Riemann Zeta function::

        sage: complex_plot(zeta, (-30,30), (-30,30))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(zeta, (-30,30), (-30,30)))

    Extra options will get passed on to show(), as long as they are valid::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3), figsize=[1,1])
        Graphics object consisting of 1 graphics primitive

    ::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3)).show(figsize=[1,1]) # These are equivalent

    REFERENCES:

    Plotting complex functions with colormaps follows the strategy
    from [LD2021]_ and incorporates contour techniques described
    in [WegSem2010]_.

    TESTS:

    Test to make sure that using fast_callable functions works::

        sage: f(x) = x^2
        sage: g = fast_callable(f, domain=CC, vars='x')
        sage: h = fast_callable(f, domain=CDF, vars='x')
        sage: P = complex_plot(f, (-10, 10), (-10, 10))
        sage: Q = complex_plot(g, (-10, 10), (-10, 10))
        sage: R = complex_plot(h, (-10, 10), (-10, 10))
        sage: S = complex_plot(exp(x)-sin(x), (-10, 10), (-10, 10))
        sage: P; Q; R; S
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive

    Test to make sure symbolic functions still work without declaring
    a variable.  (We don't do this in practice because it doesn't use
    fast_callable, so it is much slower.)

    ::

        sage: complex_plot(sqrt, (-5, 5), (-5, 5))
        Graphics object consisting of 1 graphics primitive
    """
    import matplotlib.cm
    import numpy as np
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    from sage.ext.fast_callable import fast_callable
    from sage.rings.complex_double import CDF

    try:
        f = fast_callable(f, domain=CDF, expect_one_var=True)
    except (AttributeError, TypeError, ValueError):
        pass

    cdef double x, y
    _, ranges = setup_for_eval_on_grid([], [x_range, y_range],
                                       options['plot_points'])
    x_range = ranges[0]
    y_range = ranges[1]
    cdef list z_values = []
    cdef list row
    for y in srange(*y_range, include_endpoint=True):
        row = []
        for x in srange(*x_range, include_endpoint=True):
            sig_check()
            row.append(f(new_CDF_element(x, y)))
        z_values.append(row)

    nrows = len(y_range)
    ncols = len(x_range)
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] rgbs = np.empty(dtype=float, shape=(nrows, ncols, 3))

    if cmap is None:
        # produce colors using the established default method
        rgbs = direct_complex_to_rgb(z_values, contoured=contoured, tiled=tiled)
    else:
        # choose colors from colormap
        if isinstance(cmap, str):
            if cmap == 'default':
                cmap = matplotlib.cm.get_cmap('turbo')
            else:
                cmap = matplotlib.cm.get_cmap(cmap)
        rgbs = cmap_complex_to_rgb(z_values, cmap=cmap, contoured=contoured, tiled=tiled)

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore=['xmin', 'xmax']))
    g.add_primitive(ComplexPlot(rgbs, x_range[:2], y_range[:2], options))
    return g
