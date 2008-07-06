r"""
Multifractal Random Walk

This module implements the fractal approach to understanding financial
markets that was pionered by Mandelbrot.  In particular, it implements
the multifractal random walk model of asset returns as developed by
Bacry, Kozhemyak, and Muzy, 2006, 'Continuous cascade models for asset
returns' and many other papers by Bacry et al. See
   \url{http://www.cmap.polytechnique.fr/~bacry/ftpPapers.html}

See also Mandelbrot's 'The Misbehavior of Markets' for a motivated
introduction to the general idea of using a self-similar approach to
modeling asset returns.

One of the main goals of this implementation is that everything by
highly optimized and ready for real world high performance simulation
work.

AUTHOR:
     -- William Stein (2008)
"""

from sage.rings.all import RDF, CDF, Integer
from sage.modules.all import vector
I = CDF.gen()

from time_series cimport TimeSeries

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double pow(double, double)
    double sqrt(double)

##################################################################
# Simultation
##################################################################

def stationary_gaussian_simulation(s, N, n=1):
    """
    Implementation of the Davies-Harte algorithm which given an
    autocovariance sequence (ACVS) s and an integer N, simulates N
    steps of the corresponding stationary Gaussian process with mean
    0. We assume that a certain Fourier transform associated to s is
    nonnegative; if it isn't, this algorithm fails with a
    NotImplementedError.

    INPUT:
        s -- a list of real numbers that defines the ACVS.
             Optimally s should have length N+1; if not
             we pad it with extra 0's until it has length N+1.
        N -- a positive integer

    OUTPUT:
        a list of n time series

    EXAMPLES:
    We define an autocovariance sequence:
        sage: N = 2^15
        sage: s = [1/math.sqrt(k+1) for k in [0..N]]
        sage: s[:5]
        [1.0, 0.70710678118654746, 0.57735026918962584, 0.5, 0.44721359549995793]

    We run the simulation:
        sage: set_random_seed(0)
        sage: sim = finance.stationary_gaussian_simulation(s, N)[0]

    Note that indeed the autocovariance sequence approximates s well:
        sage: [sim.autocovariance(i) for i in [0..4]]
        [0.98665816086255..., 0.69201577095377..., 0.56234006792017..., 0.48647965409871..., 0.43667043322102...]

    WARNING: If you were to do the above computation with a small
    value of N, then the autocovariance sequence would not approximate
    s very well.

    REFERENCES:
    This is a standard algorithm that is described in several papers.
    It is summarized nicely with many applications at the beginning of
    'SIMULATING A CLASS OF STATIONARY GAUSSIAN PROCESSES USING THE
    DAVIES-HARTE ALGORITHM, WITH APPLICATION TO LONG MEMORY
    PROCESSES', 2000, Peter F. Craigmile, which is easily found as a
    free pdf via a Google search.  This paper also generalizes the
    algorithm to the case when all elements of s are nonpositive.

    The book 'Wavelet Methods for Time Series Analysis' by Percival
    and Walden also describes this algorithm, but has a typo in that
    they put a 2*pi instead of pi a certain sum.  That book describes
    exactly how to use Fourier transform.  The description is in
    Section 7.8.  Note that these pages are missing from the Google
    Books version of the book, but are in the Amazon.com preview of
    the book.
    """
    N = Integer(N)
    if N < 1:
        raise ValueError, "N must be positive"

    if not isinstance(s, TimeSeries):
        s = TimeSeries(s)

    # Make sure s has length N+1.
    if len(s) > N + 1:
        s = s[:N+1]
    elif len(s) < N + 1:
        s += TimeSeries(N+1-len(s))

    # Make symmetrized vector.
    v = s + s[1:-1].reversed()

    # Compute its fast fourier transform.
    cdef TimeSeries a
    a = v.fft()

    # Take the real entries in the result.
    a = TimeSeries([a[0]]) + a[1:].scale_time(2)

    # Verify the nonnegativity condition.
    if a.min() < 0:
        raise NotImplementedError, "Stationary Gaussian simulation only implemented when Fourier transform is nonnegative"

    sims = []
    cdef Py_ssize_t i, k, iN = N
    cdef TimeSeries y, Z
    cdef double temp, nn = N
    for i from 0 <= i < n:
        Z = TimeSeries(2*iN).randomize('normal')
        y = TimeSeries(2*iN)
        y._values[0] = sqrt(2*nn*a._values[0]) * Z._values[0]
        y._values[2*iN-1] = sqrt(2*nn*a._values[iN]) * Z._values[2*iN-1]
        for k from 1 <= k < iN:
            temp = sqrt(nn*a._values[k])
            y._values[2*k-1] = temp*Z._values[2*k-1]
            y._values[2*k] = temp*Z._values[2*k]
        y.ifft(overwrite=True)
        sims.append(y[:iN])

    return sims

def fractional_gaussian_noise_simulation(double H, double sigma2, N, n=1):
    """
    Return $n$ simulations of with N steps each of fractional Gaussian
    noise with Hurst parameter $H$ and innovations variance sigma2.

    INPUT:
        H -- float; 0 < H < 1; the Hurst parameter
        sigma2 - float; innovation variance
        N -- positive integer
        n -- positive integer (default: 1)

    OUTPUT:
        list of n time series.

    EXAMPLES:
    We simulate a fractional Gaussian noise.
        sage: set_random_seed(0)
        sage: finance.fractional_gaussian_noise_simulation(0.8,1,10,2)
        [[-0.1157, 0.7025, 0.4949, 0.3324, 0.7110, 0.7248, -0.4048, 0.3103, -0.3465, 0.2964],
         [-0.5981, -0.6932, 0.5947, -0.9995, -0.7726, -0.9070, -1.3538, -1.2221, -0.0290, 1.0077]]

    The sums define a fractional Brownian motion process:
        sage: set_random_seed(0)
        sage: finance.fractional_gaussian_noise_simulation(0.8,1,10,1)[0].sums()
        [-0.1157, 0.5868, 1.0818, 1.4142, 2.1252, 2.8500, 2.4452, 2.7555, 2.4090, 2.7054]

    ALGORITHM: See 'SIMULATING A CLASS OF STATIONARY GAUSSIAN
    PROCESSES USING THE DAVIES-HARTE ALGORITHM, WITH APPLICATION TO
    LONG MEMORY PROCESSES', 2000, Peter F. Craigmile for a discussion
    and references for why the algorithm we give -- which uses
    the stationary_gaussian_simulation
    """
    if H <= 0 or H >= 1:
        raise ValueError, "H must satisfy 0 < H < 1"
    if sigma2 <= 0:
        raise ValueError, "sigma2 must be positive"
    N = Integer(N)
    if N < 1:
        raise ValueError, "N must be positive"
    cdef TimeSeries s = TimeSeries(N+1)
    s._values[0] = sigma2
    cdef Py_ssize_t k
    cdef double H2 = 2*H
    for k from 1 <= k <= N:
        s._values[k] = sigma2/2 * (pow(k+1,H2) - 2*pow(k,H2) + pow(k-1,H2))
    return stationary_gaussian_simulation(s, N, n)


def fractional_brownian_motion_simulation(double H, double sigma2, N, n=1):
    """
    Returns the partial sums of a fractional Gaussian noise simulation
    with the same input parameters.

    INPUT:
        H -- float; 0 < H < 1; the Hurst parameter
        sigma2 - float; innovation variance
        N -- positive integer
        n -- positive integer (default: 1)
    OUTPUT:
        list of n time series.

    EXAMPLES:
        sage: set_random_seed(0)
        sage: finance.fractional_brownian_motion_simulation(0.8,0.1,8,1)
        [[-0.0754, 0.2629, 0.0861, 0.2324, 0.1765, -0.0557, 0.0198, -0.0177]]
        sage: set_random_seed(0)
        sage: finance.fractional_brownian_motion_simulation(0.8,0.01,8,1)
        [[-0.0239, 0.0831, 0.0272, 0.0735, 0.0558, -0.0176, 0.0063, -0.0056]]
        sage: finance.fractional_brownian_motion_simulation(0.8,0.01,8,2)
        [[-0.0167, 0.0510, -0.0081, 0.0595, 0.0878, 0.0806, -0.1131, 0.0283],
         [0.0244, -0.0397, 0.0278, -0.0489, 0.1127, 0.0245, 0.0589, 0.0535]]
    """
    return fractional_gaussian_noise_simulation(H,sigma2,N,n)

def multifractal_cascade_random_walk_simulation(double T,
                                                double lambda2,
                                                double ell,
                                                double sigma2,
                                                N,
                                                n=1):
    r"""
    Return a list of n simulations of a multifractal random walk using
    the log-normal cascade model of Bacry-Kozhemyak-Muzy 2008.  This
    walk can be interpreted as the sequence of logarithms of a price
    series.

    INPUT:
        T       -- positive real; the integral scale
        lambda2 -- positive real; the intermittency coefficient
        ell     -- a small number -- time step size.
        sigma2  -- variance of the Gaussian white noise eps[n]
        N       -- number of steps in each simulation
        n       -- the number of separate simulations to run

    OUTPUT:
        list of time series

    EXAMPLES:
        sage: set_random_seed(0)
        sage: a = finance.multifractal_cascade_random_walk_simulation(3770,0.02,0.01,0.01,10,3)
        sage: a
        [[-0.0096, 0.0025, 0.0066, 0.0016, 0.0078, 0.0051, 0.0047, -0.0013, 0.0003, -0.0043],
         [0.0003, 0.0035, 0.0257, 0.0358, 0.0377, 0.0563, 0.0661, 0.0746, 0.0749, 0.0689],
         [-0.0120, -0.0116, 0.0043, 0.0078, 0.0115, 0.0018, 0.0085, 0.0005, 0.0012, 0.0060]]

    The corresponding price series:
        sage: a[0].exp()
        [0.9905, 1.0025, 1.0067, 1.0016, 1.0078, 1.0051, 1.0047, 0.9987, 1.0003, 0.9957]

    MORE DETAILS: The random walk has $n$th step $\eps_n
    e^{\omega_n}$, where $\eps_n$ is gaussian white noise of variance
    $\sigma^2$ and $\omega_n$ is renormalized gaussian magnitude,
    which is given by a stationary gaussian simultation associated to
    a certain autocovariance sequence.  See Bacry, Kozhemyak, Muzy,
    2006, 'Continuous cascade models for asset returns' for details.

    """
    if ell <= 0:
        raise ValueError, "ell must be positive"
    if T <= ell:
        raise ValueError, "T must be > ell"
    if lambda2 <= 0:
        raise ValueError, "lambda2 must be positive"
    N = Integer(N)
    if N < 1:
        raise ValueError, "N must be positive"

    # Compute the mean of the Gaussian stationary process omega.
    # See page 3 of Bacry, Kozhemyak, Muzy, 2008 -- "Log-Normal
    # Continuous Cascades..."
    cdef double mean = -lambda2 * (log(T/ell) + 1)

    # Compute the autocovariance sequence of the Gaussian
    # stationary process omega.  [Loc. cit.]  We have
    # tau = ell*i in that formula (6).
    cdef TimeSeries s = TimeSeries(N+1)
    s._values[0] = lambda2*(log(T/ell) + 1)
    cdef Py_ssize_t i
    for i from 1 <= i < N+1:
        if ell*i >= T:
            s._values[i] = lambda2 * log(T/(ell*i))
        else:
            break  # all covariance numbers after this point are 0

    # Compute n simultations of omega, but with mean 0.
    omega = stationary_gaussian_simulation(s, N+1, n)

    # Increase each by the given mean.
    omega = [t.add_scalar(mean) for t in omega]

    # For each simultation, create corresponding multifractal
    # random walk, as explained on page 6 of [loc. cit.]
    sims = []
    cdef Py_ssize_t k
    cdef TimeSeries eps, steps, om
    for om in omega:
        # First compute N Gaussian white noise steps
        eps = TimeSeries(N).randomize('normal', 0, sigma2)

        # Compute the steps of the multifractal random walk.
        steps = TimeSeries(N)
        for k from 0 <= k < N:
            steps._values[k] = eps._values[k] * exp(om._values[k])

        sims.append(steps.sums())

    return sims







## ##################################################################
## # Forecasting
## ##################################################################
## def multifractal_cascade_linear_filter(double T,  double lambda2,
##                    double sigma, double ell, Py_ssize_t M):
##     """
##     NOTE: I CANT GET ANALYTIC PARAMETERS RIGHT -- SWITCH TO MONTE CARLO...

##     Create the linear filter for the multifractal cascade with M steps
##     with given multifractal parameters.  Use this to predict the
##     squares of log price returns, i.e., the volatility.

##     INPUT:
##         T       -- positive real; the integral scale
##         lambda2 -- positive real; the intermittency coefficient
##         sigma   -- scaling parameter
##         ell     -- a small number -- time step size.
##         M       -- number of steps in linear filter

##     OUTPUT:
##         a time series, which should be given as the input
##         to the linear_forecast function.

##     REFERENCES:  See Section 13.6 "Linear Prediction and Linear Predictive Coding"
##     from Numerical Recipes in C for computing linear prediction
##     filters.  See the top of page 3 -- equation 39 -- of Bacry,
##     Kozhemyak, Muzy, 2006, Continuous cascade models for asset returns
##     for the specific analytic formulas we use of the theoretical
##     autocovariance.

##     ALGORITHM: The runtime is dominated by solving a numerical linear
##     equation involving an M x M matrix. Thus the filter for M up to
##     5000 can reasonably be done in less than a minute.

##     EXAMPLES:
##     We compute a linear filter with just a few terms:
##         sage: L = finance.multifractal_cascade_linear_filter(1000, 0.02, 2, 0.01, 5); L
##         [0.6453, 0.0628, 0.0822, 0.0570, 0.0959]

##     Note that the sum as we increase the length of the filter goes to 1:
##         sage: sum(L)
##         0.94318436639725745
##         sage: L = finance.multifractal_cascade_linear_filter(1000, 0.02, 2, 0.01, 1000); L
##         [0.6014, 0.0408, 0.0576, 0.0323, 0.0241 ... 0.0001, 0.0001, 0.0002, 0.0002, 0.0005]
##         sage: sum(L)
##         0.99500051396305267

##     We create a random multifractal walk:
##         sage: set_random_seed(0)
##         sage: y = finance.multifractal_cascade_random_walk_simulation(3770,0.02,0.01,2000,1)[0]; y
##         [0.0002, -0.0035, -0.0107, 0.0234, 0.0181 ... 0.2526, 0.2427, 0.2508, 0.2530, 0.2530]

##     We then square each term in the series.
##         sage: y2 = y.pow(2)
##         sage: y2
##         [0.0000, 0.0000, 0.0001, 0.0005, 0.0003 ... 0.0638, 0.0589, 0.0629, 0.0640, 0.0640]
##         sage: y2[:-1].linear_forecast(L)
##         0.064481359234932076
##         sage: y2[:-2].linear_forecast(L)
##         0.063950875470110996

##     """
##     cdef TimeSeries c
##     cdef Py_ssize_t i

##     if M <= 0:
##         raise ValueError, "M must be positive"

##     # 1. Compute the autocovariance sequence
##     c = TimeSeries(M+1)

##     for i from 1 <= i <= M:
##         c._values[i] = lambda2 * log(T/i)
##     print c

##     # cdef double s4 = pow(sigma,4)
##     # cdef double ell2 = ell*ell
##     # cdef e = -4*lambda2
##     #for i from 1 <= i <= M:
##     #   c(i) = sigma^4 * ell^2 * (i/T)^(-4*lambda^2)
##     #    c._values[i] = s4 * ell2 * pow(i/T, e)

##     #cdef double K, tau
##     #for i from 0 <= i <= M:
##     #    tau = ell*i
##     #    K = s4 * pow(T,4*lambda2) / ((1 + e) * (2 + e))
##     #    c._values[i] = K * (pow(abs(ell + tau), 2+e) + pow(abs(ell-tau), 2+e) - 2*pow(abs(tau),2+e))

##     # Also create the numpy vector v with entries c(1), ..., c(M).
##     v = c[1:].numpy()

##     # 2. Create the autocovariances numpy 2d array A whose i,j entry
##     # is c(|i-j|).
##     import numpy
##     A = numpy.array([[c[abs(j-k)] for j in range(M)] for k in range(M)])

##     # 3. Solve the equation A * x = v
##     x = numpy.linalg.solve(A, v)

##     # 4. The entries of x give the linear filter coefficients.
##     return TimeSeries(x)



##################################################################
# Parameter Estimation
##################################################################
