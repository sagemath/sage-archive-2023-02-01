"""
Riemann Mapping

AUTHORS:

- Ethan Van Andel (2009): initial version

- Robert Bradshaw (2009): his "complex_plot" was adapted for plot_colored

Development supported by NSF award No. 0702939.
"""

#*****************************************************************************
#       Copyright (C) 2009 Ethan Van Andel <evlutte@gmail.com>,
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

include "../ext/stdsage.pxi"
include "../ext/interrupt.pxi"

from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options
from sage.plot.all import list_plot, Graphics
from sage.plot.misc import setup_for_eval_on_grid

from sage.ext.fast_eval import fast_callable

from sage.rings.all import CDF

from sage.misc.misc import srange

from sage.gsl.interpolation import spline

import numpy as np
cimport numpy as np

from math import pi
from math import sin
from math import cos

from cmath import exp
from cmath import phase

cdef double PI = pi
cdef double TWOPI = 2*PI
cdef I = complex(0,1)

FLOAT = np.float64
ctypedef np.float64_t FLOAT_T


cdef class Riemann_Map:
    """
    The ``Riemann_Map`` class computes a Riemann or Ahlfors map from
    supplied data. It also has various methods to provide information about
    the map. A Riemann map conformally maps a simply connected region in
    the complex plane to the unit disc. The Ahlfors map does the same thing
    for multiply connected regions.

    Note that all the methods are numeric rather than analytic, for unusual
    regions or insufficient collocation points may give very inaccurate
    results.

    INPUT:

    - ``fs`` -- A list of the boundaries of the region, given as
      complex-valued functions with domain ``0`` to ``2*pi``. Note that the
      outer boundary must be parameterized counter clockwise
      (i.e. ``e^(I*t)``) while the inner boundaries must be clockwise
      (i.e. ``e^(-I*t)``).

    - ``fprimes`` -- A list of the derivatives of the boundary functions.
      Must be in the same order as ``fs``.

    - ``a`` -- Complex, the center of the Riemann map. Will be mapped to
      the origin of the unit disc.

    The following inputs must all be passed in as named parameters:

    - ``N`` -- integer (default: ``500``), the number of collocation points
      used to compute the map. More points will give more accurate results,
      especially near the boundaries, but will take longer to compute.

    - ``ncorners`` -- integer (default: ``4``), if mapping a figure with
      (equally t-spaced) corners, better results may be obtained by
      accurately giving this parameter. Used to add the proper constant to
      the theta correspondance function.

    - ``opp`` -- boolean (default: ``False``), set to ``True`` in very rare
      cases where the theta correspondance function is off by ``pi``, that
      is, if red is mapped left of the origin in the color plot.

    EXAMPLES:

    The unit circle identity map::

        sage: m = Riemann_Map([lambda t: e^(I*t)], [lambda t: I*e^(I*t)], 0)  # long time (4 sec)

    The unit circle with a small hole::

        sage: f(t) = e^(I*t)
        sage: fprime(t) = I*e^(I*t)
        sage: hf(t) = 0.5*e^(-I*t)
        sage: hfprime(t) = 0.5*-I*e^(-I*t)
        sage: m = Riemann_Map([f, hf], [hf, hfprime], 0.5 + 0.5*I)

    A square::

        sage: ps = polygon_spline([(-1, -1), (1, -1), (1, 1), (-1, 1)])  # long time
        sage: f = lambda t: ps.value(real(t))  # long time
        sage: fprime = lambda t: ps.derivative(real(t))  # long time
        sage: m = Riemann_Map([f], [fprime], 0.25, ncorners=4)  # long time
        sage: m.plot_colored() + m.plot_spiderweb()  # long time

    Compute rough error for this map::

        sage: x = 0.75  # long time
        sage: print "error =", m.inverse_riemann_map(m.riemann_map(x)) - x  # long time
        error = (-0.000...+0.0016...j)

    ALGORITHM:

    This class computes the Riemann Map via the Szego kernel using an
    adaptation of the method described by [KT]_.

    REFERENCES:

    .. [KT] N. Kerzman and M. R. Trummer. "Numerical Conformal Mapping via
      the Szego kernel". Journal of Computational and Applied Mathematics,
      14(1-2): 111--123, 1986.
    """
    cdef int N, B, ncorners
    cdef f
    cdef opp
    cdef double complex a
    cdef np.ndarray tk, tk2, cps, dps, szego, p_vector, pre_q_vector
    cdef np.ndarray p_vector_inverse, sinalpha, cosalpha, theta_array
    cdef x_range, y_range

    def __init__(self, fs, fprimes, a, int N=500, int ncorners=4, opp=False):
        """
        Initializes the ``Riemann_Map`` class. See the class ``Riemann_Map``
        for full documentation on the input of this initialization method.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
        """
        a = np.complex128(a)
        if N <= 2:
            raise ValueError(
                "The number of collocation points must be > 2.")
        try:
            fs = [fast_callable(f, domain=CDF) for f in fs]
            fprimes = [fast_callable(f, domain=CDF) for f in fprimes]
        except AttributeError:
            pass
        self.f = fs[0]
        self.a = a
        self.ncorners = ncorners
        self.N = N  # Number of collocation pts
        self.opp = opp
        cdef int i, k
        self.tk = np.array(np.arange(N) * TWOPI / N + 0.001 / N,
                           dtype=np.float64)
        self.tk2 = np.zeros(N + 1, dtype=np.float64)
        for i in xrange(N):
            self.tk2[i] = self.tk[i]
        self.tk2[N] = TWOPI
        self.B = len(fs) # number of boundaries of the figure
        self.cps = np.zeros([self.B, N], dtype=np.complex128)
        self.dps = np.zeros([self.B, N], dtype=np.complex128)
        # Find the points on the boundaries and their derivatives.
        for k in xrange(self.B):
            for i in xrange(N):
                self.cps[k, i] = np.complex(fs[k](self.tk[i]))
                self.dps[k, i] = np.complex(fprimes[k](self.tk[i]))
        cdef double xmax = self.cps.real.max()
        cdef double xmin = self.cps.real.min()
        cdef double ymax = self.cps.imag.max()
        cdef double ymin = self.cps.imag.min()
        cdef double space = 0.1 * max(xmax - xmin, ymax - ymin)
        #The default plotting window, mainly for color plot.
        self.x_range = (xmin - space, xmax + space)
        self.y_range = (ymin - space, ymax + space)
        self._generate_theta_array()
        self._generate_interior_mapper()
        self._generate_inverse_mapper()

    def _repr_(self):
        """
        Return a string representation of this ``Riemann_Map`` object.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: isinstance(Riemann_Map([f], [fprime], 0)._repr_(), str)  # long time
            True
        """
        return "A Riemann mapping of a figure to the unit circle."

    cdef _generate_theta_array(self):
        """
        Generates the essential data for the Riemann map, primarily the
        Szego kernel and boundary correspondance.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
        """
        cp = self.cps.flatten()
        dp = self.dps.flatten()
        cdef int N = self.N
        cdef int NB = N * self.B
        cdef int B = self.B
        cdef int i, k
        cdef FLOAT_T saa, t0
        cdef np.ndarray[FLOAT_T, ndim=1] adp
        cdef np.ndarray[FLOAT_T, ndim=1] sadp
        cdef np.ndarray h
        cdef np.ndarray hconj
        cdef np.ndarray g
        cdef np.ndarray K
        cdef np.ndarray phi
        cdef np.ndarray[FLOAT_T, ndim=2] theta_array
        # Seting things up to use the Nystrom method
        adp = abs(dp)
        sadp = np.sqrt(adp)
        h = 1 / (TWOPI * I) * ((dp / adp) / (self.a - cp))
        hconj = np.array(map(np.complex.conjugate, h), dtype=np.complex128)
        g = -sadp * hconj
        normalized_dp=dp/adp
        C = I / N * sadp # equivalent to -TWOPI / N * 1 / (TWOPI * I) * sadp
        olderr = np.geterr()['invalid'] # checks the current error handling
        np.seterr(invalid='ignore')
        K = np.array(
            [C * sadp[t] *
             (normalized_dp/(cp-cp[t]) - (normalized_dp[t]/(cp-cp[t])).conjugate())
              for t in np.arange(NB)], dtype=np.complex128)
        np.seterr(invalid=olderr) # resets the error handling
        for i in xrange(NB):
            K[i, i] = 1
        phi = np.linalg.solve(K, g) / NB * TWOPI  # Nystrom
        # the all-important Szego kernel
        szego = np.array(phi.flatten() / np.sqrt(dp), dtype=np.complex128)
        self.szego = szego.reshape([B, N])
        start = 0
        # Finding the theta correspondance using phase. Misbehaves for some
        # regions.
        if B != 1:
            theta_array = np.zeros([1, NB])
            for i in xrange(NB):
                theta_array[0, i] = phase(-I * np.power(phi[i], 2) * dp[i])
            self.theta_array = np.concatenate(
                [theta_array.reshape([B, N]), np.zeros([B, 1])], axis=1)
            for k in xrange(B):
                self.theta_array[k, N] = self.theta_array[k, 0] + TWOPI
        # Finding the theta correspondance using ab. Well behaved, but
        # doesn't work on multiply connected domains.
        else:
            phi2 = phi.reshape([self.B, N])
            theta_array = np.zeros([B, N + 1], dtype=np.float64)
            for k in xrange(B):
                phik = phi2[k]
                saa = (np.dot(abs(phi), abs(phi))) * TWOPI / NB
                theta_array[k, 0] = 0
                for i in xrange(1, N):
                    theta_array[k, i] = (
                        theta_array[k, i - 1] +
                        ((TWOPI / NB * TWOPI *
                          abs(np.power(phi[1 * i], 2)) / saa +
                          TWOPI / NB * TWOPI *
                          abs(np.power(phi[1 * i - 1], 2)) / saa)) / 2)
                tmax = int(0.5 * N / self.ncorners)
                # Finding the initial value of the theta function.
                phimax = -I * phik[tmax]**2 * self.dps[k, tmax]
                if self.opp:
                    t0 = theta_array[k, tmax] + phase(phimax)
                else:
                    t0 = theta_array[k, tmax] - phase(phimax)
                for i in xrange(N):
                    theta_array[k, i] = theta_array[k, i] - t0
                theta_array[k, N] = TWOPI + theta_array[k, 0]
            self.theta_array = theta_array

    def get_szego(self, int boundary=-1, absolute_value=False):
        """
        Returns a discretized version of the Szego kernel for each boundary
        function.

        INPUT:

        The following inputs must all be passed in as named parameters:

        - ``boundary`` -- integer (default: ``-1``) if < 0,
          ``get_theta_points()`` will return the points for all boundaries.
          If >= 0, ``get_theta_points()`` will return only the points for
          the boundary specified.

        - ``absolute_value`` -- boolean (default: ``False``) if ``True``, will
          return the absolute value of the (complex valued) Szego kernel
          instead of the kernel itself. Useful for plotting.

        OUTPUT:

        A list of points of the form
        ``[t value, value of the Szego kernel at that t]``.

        EXAMPLES:

        Generic use::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: sz = m.get_szego(boundary=0)
            sage: points = m.get_szego(absolute_value=True)
            sage: list_plot(points)

        Extending the points by a spline::

            sage: s = spline(points)
            sage: s(3*pi / 4)
            0.00121587378429...

        The unit circle with a small hole::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: hf(t) = 0.5*e^(-I*t)
            sage: hfprime(t) = 0.5*-I*e^(-I*t)
            sage: m = Riemann_Map([f, hf], [hf, hfprime], 0.5 + 0.5*I)

        Getting the szego for a specifc boundary::

            sage: sz0 = m.get_szego(boundary=0)
            sage: sz1 = m.get_szego(boundary=1)
        """
        cdef int k, B
        if boundary < 0:
            temptk = self.tk
            for i in xrange(self.B - 1):
                temptk = np.concatenate([temptk, self.tk])
            if absolute_value:
                return np.column_stack(
                    [temptk, abs(self.szego.flatten())]).tolist()
            else:
                return np.column_stack([temptk, self.szego.flatten()]).tolist()
        else:
            if absolute_value:
                return np.column_stack(
                    [self.tk, abs(self.szego[boundary])]).tolist()
            else:
                return np.column_stack(
                    [self.tk, self.szego[boundary]]).tolist()

    def get_theta_points(self, int boundary=-1):
        """
        Returns an array of points of the form
        ``[t value, theta in e^(I*theta)]``, that is, a discretized version
        of the theta/boundary correspondence function. For multiply
        connected domains, ``get_theta_points`` will list the points for
        each boundary in the order that they were supplied.

        INPUT:

        The following input must all be passed in as named parameters:

        - ``boundary`` -- integer (default: ``-``1) if < 0,
          ``get_theta_points()`` will return the points for all boundaries.
          If >= 0, ``get_theta_points()`` will return only the points for
          the boundary specified.

        OUTPUT:

        A list of points of the form ``[t value, theta in e^(I*theta)]``.

        EXAMPLES:

        Getting the list of points, extending it via a spline, getting the
        points for only the outside of a multiply connected domain::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: points = m.get_theta_points()
            sage: list_plot(points)

        Extending the points by a spline::

            sage: s = spline(points)
            sage: s(3*pi / 4)
            1.62766037996...

        The unit circle with a small hole::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: hf(t) = 0.5*e^(-I*t)
            sage: hfprime(t) = 0.5*-I*e^(-I*t)
            sage: m = Riemann_Map([f, hf], [hf, hfprime], 0.5 + 0.5*I)

        Getting the szego for a specifc boundary::

            sage: tp0 = m.get_theta_points(boundary=0)
            sage: tp1 = m.get_theta_points(boundary=1)
        """
        if boundary < 0:
            temptk = self.tk2
            for i in xrange(self.B - 1):
                temptk = np.concatenate([temptk, self.tk2])
            return np.column_stack(
                [temptk, self.theta_array.flatten()]).tolist()
        else:
            return np.column_stack(
                [self.tk2, self.theta_array[boundary]]).tolist()

    cdef _generate_interior_mapper(self):
        """
        Generates the data necessary to use the ``reimann_map()`` function.
        As much setup as possible is done here to minimize what has to be
        done in ``riemann_map()``.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
        """
        cdef int N = self.N
        cdef double complex coeff = 3*I / (8*N)
        dps = self.dps
        theta_array = self.theta_array
        cdef np.ndarray[double complex, ndim=2] p_vector = np.zeros(
            [self.B, N + 1], dtype=np.complex128)
        cdef int k, i
        # Lots of setup for Simpson's method of integration.
        for k in xrange(self.B):
            for i in xrange(N // 3):
                p_vector[k, 3*i] = (2*coeff * dps[k, 3*i] *
                                    exp(I * theta_array[k, 3*i]))
                p_vector[k, 3*i + 1] = (3*coeff * dps[k, 3*i + 1] *
                                        exp(I * theta_array[k, 3*i + 1]))
                p_vector[k, 3*i + 2] = (3*coeff * dps[k, 3*i + 2] *
                                        exp(I * theta_array[k, 3*i + 2]))
            p_vector[k, 0] = 1*coeff * dps[k, 0] * exp(I * theta_array[k, 0])
            if N % 3 == 0:
                p_vector[k, N] = 1*coeff * dps[k, 0] * exp(I*theta_array[k, 0])
            elif (N - 2) % 3 == 0:
                p_vector[k, N - 2] = ((coeff + I/(3*N)) * dps[k, N- 2 ] *
                                      exp(I * theta_array[k, N - 2]))
                p_vector[k, N- 1 ] = (4*I / (3*N) * dps[k, N - 1] *
                                      exp(I * theta_array[k, N - 1]))
                p_vector[k, N] = (I / (3*N) * dps[k, 0] *
                                  exp(I * theta_array[k, 0]))
            else:
                p_vector[k, N - 4] = ((coeff + I / (3*N)) * dps[k, N - 4] *
                                      exp(I * theta_array[k, N - 4]))
                p_vector[k, N - 3] = (4*I / (3*N) * dps[k, N - 3] *
                                      exp(I * theta_array[k, N - 3]))
                p_vector[k, N - 2] = (2*I / (3*N) * dps[k, N - 2] *
                                      exp(I * theta_array[k, N - 2]))
                p_vector[k, N - 1] = (4*I / (3*N) * dps[k, N - 1] *
                                      exp(I * theta_array[k, N - 1]))
                p_vector[k, N] = (I / (3*N) * dps[k, 0] *
                                  exp(I * theta_array[k, 0]))
        self.p_vector = p_vector.flatten()
        cdef np.ndarray[double complex, ndim=1] pq = self.cps[:,list(range(N))+[0]].flatten()
        self.pre_q_vector = pq

    cpdef riemann_map(self, pt):
        """
        Returns the Riemann mapping of a point. That is, given ``pt`` on
        the interior of the mapped region, ``riemann_map()`` will return
        the point on the unit disk that ``pt`` maps to. Note that this
        method only works for interior points; it breaks down very close
        to the boundary. To get boundary corrospondance, use
        ``get_theta_points()``.

        INPUT:

        - ``pt`` -- A complex number representing the point to be
          inverse mapped.

        OUTPUT:

        A complex number representing the point on the unit circle that
        the input point maps to.

        EXAMPLES:

        Can work for different types of complex numbers::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: m.riemann_map(0.25 + sqrt(-0.5))
            (0.137514...+0.87669602...j)
            sage: m.riemann_map(1.3*I)
            (-1.56...e-05+0.989694...j)
            sage: I = CDF.gen()
            sage: m.riemann_map(0.4)
            (0.733242677...+3.2...e-06j)
            sage: import numpy as np
            sage: m.riemann_map(np.complex(-3, 0.0001))
            (1.405757...e-05+8.06...e-10j)
        """
        pt1 = np.complex(pt)
        cdef np.ndarray[double complex, ndim=1] q_vector = 1 / (
            self.pre_q_vector - pt1)
        return -np.dot(self.p_vector, q_vector)

    cdef _generate_inverse_mapper(self):
        """
        Generates the data necessary to use the
        ``inverse_reimann_map()`` function. As much setup as possible is
        done here to minimize what has to be done in
        ``inverse_riemann_map()``.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
        """
        cdef int N = self.N
        cdef int B = self.B
        cdef float di
        theta_array = self.theta_array
        self.p_vector_inverse = np.zeros([B, N], dtype=np.complex128)
        # Setup for trapezoid integration because integration points are
        # not equally spaced.
        for k in xrange(B):
            for i in xrange(N):
                di = theta_array[k, (i + 1) % N] - theta_array[k, (i - 1) % N]
                if di > PI:
                    di = di - TWOPI
                elif di < -PI:
                    di = di + TWOPI
                self.p_vector_inverse[k, i] = di / 2
        self.sinalpha = np.zeros([B, N], dtype=np.float64)
        for k in xrange(B):
            for i in xrange(N):
                self.sinalpha[k, i] = sin(-theta_array[k, i])
        self.cosalpha = np.zeros([B, N], dtype=np.float64)
        for k in xrange(B):
            for i in xrange(N):
                self.cosalpha[k, i] = cos(-theta_array[k, i])

    def inverse_riemann_map(self, pt):
        """
        Returns the inverse Riemann mapping of a point. That is, given ``pt``
        on the interior of the unit disc, ``inverse_reimann_map()`` will
        return the point on the original region that would be Riemann
        mapped to ``pt``.

        .. NOTE::

            This method does not work for multiply connected domains.

        INPUT:

        - ``pt`` -- A complex number (usually with absolute value <= 1)
          representing the point to be inverse mapped.

        OUTPUT:

        The point on the region that Riemann maps to the input point.

        EXAMPLES:

        Can work for different types of complex numbers::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: m.inverse_riemann_map(0.5 + sqrt(-0.5))
            (0.406880548363...+0.361470279816...j)
            sage: m.inverse_riemann_map(0.95)
            (0.486319431795...-4.90019052...j)
            sage: m.inverse_riemann_map(0.25 - 0.3*I)
            (0.165324498558...-0.180936785500...j)
            sage: import numpy as np
            sage: m.inverse_riemann_map(np.complex(-0.2, 0.5))
            (-0.156280570579...+0.321819151891...j)
        """
        pt = np.complex128(pt)
        r = abs(pt)
        if r == 0:
            stheta = 0
            ctheta = 0
        else:
            stheta = pt.imag / r
            ctheta = pt.real / r
        k = 0
        return (1 - r**2) / TWOPI * np.dot(
                self.p_vector_inverse[k] * self.cps[k],
                1 / (1 + r**2 - 2*r *
                     (ctheta * self.cosalpha[k] - stheta * self.sinalpha[k])))

    def plot_boundaries(self, plotjoined=True, rgbcolor=[0,0,0], thickness=1):
        """
        Plots the boundaries of the region for the Riemann map. Note that
        this method DOES work for multiply connected domains.

        INPUT:

        The following inputs must all be passed in as named parameters:

        - ``plotjoined`` -- boolean (default: ``True``) If ``False``,
          discrete points will be drawn; otherwise they will be connected
          by lines. In this case, if ``plotjoined=False``, the points shown
          will be the original collocation points used to generate the
          Riemann map.

        - ``rgbcolor`` -- float array (default: ``[0,0,0]``) the
          red-green-blue color of the boundary.

        - ``thickness`` -- positive float (default: ``1``) the thickness of
          the lines or points in the boundary.

        EXAMPLES:

        General usage::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)

        Default plot::

            sage: m.plot_boundaries()

        Big blue collocation points::

            sage: m.plot_boundaries(plotjoined=False, rgbcolor=[0,0,1], thickness=6)
        """
        plots = range(self.B)
        for k in xrange(self.B):
            # This should be eliminated when the thickness/pointsize issue
            # is resolved later. Same for the others in plot_spiderweb().
            if plotjoined:
                plots[k] = list_plot(
                    comp_pt(self.cps[k], 1), plotjoined=True,
                    rgbcolor=rgbcolor, thickness=thickness)
            else:
                plots[k] = list_plot(
                    comp_pt(self.cps[k], 1), rgbcolor=rgbcolor,
                    pointsize=thickness)
        return sum(plots)

    def plot_spiderweb(self, spokes=16, circles=4, pts=32, linescale=0.99,
                       rgbcolor=[0,0,0], thickness=1, plotjoined=True):
        """
        Generates a traditional "spiderweb plot" of the Riemann map. Shows
        what concentric circles and radial lines map to. Note that this
        method DOES NOT work for multiply connected domains.

        INPUT:

        The following inputs must all be passed in as named parameters:

        - ``spokes`` -- integer (default: ``16``) the number of equally
          spaced radial lines to plot.

        - ``circles`` -- integer (default: ``4``) the number of equally
          spaced circles about the center to plot.

        - ``pts`` -- integer (default: ``32``) the number of points to
          plot. Each radial line is made by ``1*pts`` points, each circle
          has ``2*pts`` points. Note that high values may cause erratic
          behavior of the radial lines near the boundaries.

        - ``linescale`` -- float between 0 and 1. Shrinks the radial lines
          away from the boundary to reduce erratic behavior.

        - ``rgbcolor`` -- float array (default: ``[0,0,0]``) the
          red-green-blue color of the spiderweb.

        - ``thickness`` -- positive float (default: ``1``) the thickness of
          the lines or points in the spiderweb.

        - ``plotjoined`` -- boolean (default: ``True``) If ``False``,
          discrete points will be drawn; otherwise they will be connected
          by lines.

        EXAMPLES:

        General usage::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)

        Default plot::

            sage: m.plot_spiderweb()

        Simplified plot with many discrete points::

            sage: m.plot_spiderweb(spokes=4, circles=1, pts=400, linescale=0.95, plotjoined=False)

        Plot with thick, red lines::

            sage: m.plot_spiderweb(rgbcolor=[1,0,0], thickness=3)

        To generate the unit circle map, it's helpful to see what the
        original spiderweb looks like::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: m = Riemann_Map([f], [fprime], 0, 1000)
            sage: m.plot_spiderweb()
        """
        edge = self.plot_boundaries(plotjoined=plotjoined, rgbcolor=rgbcolor,
                                    thickness=thickness)
        circle_list = range(circles)
        theta_array = self.theta_array[0]
        s = spline(np.column_stack([self.theta_array[0], self.tk2]).tolist())
        tmax = self.theta_array[0, self.N]
        tmin = self.theta_array[0, 0]
        cdef int k, i
        for k in xrange(circles):
            temp = range(pts*2)
            for i in xrange(2*pts):
                temp[i] = self.inverse_riemann_map(
                    (k + 1) / (circles + 1.0) * exp(I*i * TWOPI / (2*pts)))
            if plotjoined:
                circle_list[k] = list_plot(
                    comp_pt(temp, 1), rgbcolor=rgbcolor, thickness=thickness,
                    plotjoined=True)
            else:
                circle_list[k] = list_plot(
                    comp_pt(temp, 1), rgbcolor=rgbcolor, pointsize=thickness)
        line_list = range(spokes)
        for k in xrange(spokes):
            temp = range(pts)
            angle = (k*1.0) / spokes * TWOPI
            if angle >= tmax:
                angle -= TWOPI
            elif angle <= tmin:
                angle += TWOPI
            for i in xrange(pts - 1):
                temp[i] = self.inverse_riemann_map(
                    (i * 1.0) / (pts * 1.0) * exp(I * angle) * linescale)
            temp[pts - 1] = np.complex(
                self.f(s(angle)) if angle <= tmax else self.f(s(angle-TWOPI)))
            if plotjoined:
                line_list[k] = list_plot(
                    comp_pt(temp, 0), rgbcolor=rgbcolor, thickness=thickness,
                    plotjoined=True)
            else:
                line_list[k] = list_plot(
                    comp_pt(temp, 0), rgbcolor=rgbcolor, pointsize=thickness)
        return edge + sum(circle_list) + sum(line_list)

    def plot_colored(self, plot_range=[], int plot_points=100):
        """
        Draws a colored plot of the Riemann map. A red point on the
        colored plot corresponds to a red point on the unit disc. Note that
        this method DOES work for multiply connected domains.

        INPUT:

        The following inputs must all be passed in as named parameters:

        - ``plot_range`` -- (default: ``[]``) list of 4 values
          ``(xmin, xmax, ymin, ymax)``. Declare if you do not want the plot
          to use the default range for the figure.

        - ``plot_points`` -- integer (default: ``100``), number of points
          to plot in each direction of the grid. Note that very large values
          can cause this function to run slowly.

        EXAMPLES:

        Given a Riemann map m, general usage::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: m.plot_colored()

        Plot zoomed in on a specific spot::

            sage: m.plot_colored(plot_range=[-1,1,2,3])

        High resolution plot::

            sage: m.plot_colored(plot_points=1000)  # long time (30s on sage.math, 2011)

        To generate the unit circle map, it's helpful to see what the
        colors correspond to::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: m = Riemann_Map([f], [fprime], 0, 1000)
            sage: m.plot_colored()
        """
        cdef int i, j
        cdef double xmin
        cdef double xmax
        cdef double ymin
        cdef double ymax
        if plot_range == []:
            xmin = self.x_range[0]
            xmax = self.x_range[1]
            ymin = self.y_range[0]
            ymax = self.y_range[1]
        else:
            xmin = plot_range[0]
            xmax = plot_range[1]
            ymin = plot_range[2]
            ymax = plot_range[3]
        xstep = (xmax - xmin) / plot_points
        ystep = (ymax - ymin) / plot_points
        cdef np.ndarray z_values = np.empty(
            [plot_points, plot_points], dtype=np.complex128)
        for i in xrange(plot_points):
            for j in xrange(plot_points):
                z_values[j, i] = self.riemann_map(
                    np.complex(xmin + 0.5*xstep + i*xstep,
                               ymin + 0.5*ystep + j*ystep))
        g = Graphics()
        g.add_primitive(ColorPlot(z_values, (xmin, xmax), (ymin, ymax)))
        return g

cdef comp_pt(clist, loop=True):
    """
    This function converts a list of complex numbers to the plottable
    `(x,y)` form. If ``loop=True``, then the first point will be
    added as the last, i.e. to plot a closed circle.

    INPUT:

    - ``clist`` -- a list of complex numbers.

    - ``loop`` -- boolean (default: ``True``) controls whether or not the
      first point will be added as the last to plot a closed circle.

    EXAMPLES:

    This tests it implicitly::

        sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
        sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
        sage: m = Riemann_Map([f], [fprime], 0)
        sage: m.plot_spiderweb()
    """
    list2 = range(len(clist) + 1) if loop else range(len(clist))
    for i in xrange(len(clist)):
        list2[i] = (clist[i].real, clist[i].imag)
    if loop:
        list2[len(clist)] = list2[0]
    return list2

cdef inline double mag_to_lightness(double r):
    """
    Tweak this to adjust how the magnitude affects the color.

    INPUT:

    - ``r`` -- a non-negative real number.

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 10]])
        array([[[ 0.        ,  0.        ,  0.        ],
                [ 0.77172568,  0.        ,  0.        ],
                [ 1.        ,  0.22134776,  0.22134776]]])
    """
    return 1 - r

cdef complex_to_rgb(np.ndarray z_values):
    r"""
    Convert from an array of complex numbers to its corresponding matrix of
    RGB values.

    INPUT:

    - ``z_values`` -- A grid of complex numbers, as a list of lists.

    OUTPUT:

    An `N \times M \times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r,g,b) tuple.

    EXAMPLES::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 1000]])
        array([[[ 0.        ,  0.        ,  0.        ],
                [ 0.77172568,  0.        ,  0.        ],
                [ 1.        ,  0.64421177,  0.64421177]]])
        sage: complex_to_rgb([[0, 1j, 1000j]])
        array([[[ 0.        ,  0.        ,  0.        ],
                [ 0.38586284,  0.77172568,  0.        ],
                [ 0.82210588,  1.        ,  0.64421177]]])
    """
    cdef unsigned int i, j, imax, jmax
    cdef double x, y, mag, arg
    cdef double lightness, hue, top, bot
    cdef double r, g, b
    cdef int ihue
    cdef z

    from cmath import phase

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef np.ndarray[np.float_t, ndim=3, mode="c"] rgb = np.empty(
        dtype=np.float, shape=(imax, jmax, 3))

    sig_on()
    for i from 0 <= i < imax:
        row = z_values[i]
        for j from 0 <= j < jmax:
            z = row[j]
            mag = abs(z)
            arg = phase(z)
            lightness = mag_to_lightness(mag)
            # in hsv, variable value, full saturation (s=1, v=1+lightness)
            if lightness < 0:
                bot = 0
                top = (1 + lightness)
            # in hsv, variable saturation, full value (v=1, s=1-lightness)
            else:
                bot = lightness
                top = 1
            # Note that does same thing as colorsys module hsv_to_rgb for
            # this setup, but in Cython.
            hue = 3*arg / PI
            # usual hsv hue is thus h=arg/(2*pi) for positive,
            # h=arg/(2*PI)+1 for negative
            if hue < 0:
                hue += 6
            ihue = <int>hue
            if ihue == 0:
                r = top
                g = bot + hue * (top - bot)
                b = bot
            elif ihue == 1:
                r = bot + (2 - hue) * (top - bot)
                g = top
                b = bot
            elif ihue == 2:
                r = bot
                g = top
                b = bot + (hue - 2) * (top - bot)
            elif ihue == 3:
                r = bot
                g = bot + (4 - hue) * (top - bot)
                b = top
            elif ihue == 4:
                r = bot + (hue - 4) * (top - bot)
                g = bot
                b = top
            else:
                r = top
                g = bot
                b = bot + (6 - hue) * (top - bot)
            rgb[i, j, 0] = r
            rgb[i, j, 1] = g
            rgb[i, j, 2] = b
    sig_off()
    return rgb

class ColorPlot(GraphicPrimitive):
    """
    The GraphicsPrimitive to display complex functions in using the domain
    coloring method

    INPUT:

        - ``z_values`` -- An array of complex values to be plotted.

        - ``x_range`` -- A minimum and maximum x value for the plot.

        - ``y_range`` -- A minimum and maximum y value for the plot.

    EXAMPLES::

        sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
    """
    def __init__(self, z_values, x_range, y_range):
        """
        Setup a ``ColorPlot`` object.

        TESTS::

            sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
        """
        self.x_range = x_range
        self.y_range = y_range
        self.z_values = z_values
        self.x_count = len(z_values)
        self.y_count = len(z_values[0])
        self.rgb_data = complex_to_rgb(z_values)
        GraphicPrimitive.__init__(self, [])

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: p = complex_plot(lambda z: z, (-1, 2), (-3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', -1.0), ('ymax', 4.0), ('ymin', -3.0)]
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.x_range, self.y_range, dict=True)

    def _allowed_options(self):
        """
        Return a dictionary of valid options for this ``ColorPlot`` object.

        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._allowed_options(), dict)
            True
        """
        return {'plot_points': 'How many points to use for plotting precision',
                'interpolation': 'What interpolation method to use'}

    def _repr_(self):
        """
        Return a string representation of this ``ColorPlot`` object.

        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._repr_(), str)
            True
        """
        return "ColorPlot defined by a %s x %s data grid" % (
            self.x_count, self.y_count)

    def _render_on_subplot(self, subplot):
        """
        Render the graphics object on a subplot.

        TESTS::

            sage: complex_plot(lambda x: x^2, (-5, 5), (-5, 5))
        """
        options = self.options()
        x0, x1 = float(self.x_range[0]), float(self.x_range[1])
        y0, y1 = float(self.y_range[0]), float(self.y_range[1])
        subplot.imshow(self.rgb_data, origin='lower',
                       extent=(x0,x1,y0,y1), interpolation='catrom')
