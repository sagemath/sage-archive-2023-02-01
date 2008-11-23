"""
Complex Elliptic Curve L-series

"""

from sage.structure.sage_object import SageObject
from sage.rings.all import (
    RealField,
    RationalField,
    ComplexField)
from math import sqrt, exp, ceil
import sage.functions.transcendental as transcendental
R = RealField()
Q = RationalField()
C = ComplexField()
import sage.functions.constants as constants
import sage.misc.all as misc

class Lseries_ell(SageObject):
    """
    An elliptic curve $L$-series.

    EXAMPLES:

    """
    def __init__(self, E):
        """
        Create an elliptic curve $L$-series.

        EXAMPLES:
            sage: EllipticCurve([1..5]).lseries()
            Complex L-series of the Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field
        """
        self.__E = E

    def elliptic_curve(self):
        """
        Return the elliptic curve that this L-series is attached to.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: L = E.lseries()
            sage: L.elliptic_curve ()
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return self.__E

    def taylor_series(self, a=1, prec=53, series_prec=6, var='z'):
        """
        Return the Taylor series of this $L$-series about $a$ to
        the given precision (in bits) and the number of terms.

        The output is a series in var, where you should view var as
        equal to s-a.  Thus this function returns the formal power
        series whose coefficients are L^{(n)}(a)/n!.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: L = E.lseries()
            sage: L.taylor_series(series_prec=3)      # random nearly 0 constant and linear terms
            -2.69129566562797e-23 + (1.52514901968783e-23)*z + 0.759316500288427*z^2 + O(z^3)
            sage: L.taylor_series(series_prec=3)[2:]
            0.759316500288427*z^2 + O(z^3)
        """
        D = self.dokchitser(prec)
        return D.taylor_series(a, series_prec, var)

    def _repr_(self):
        """
        Return string representation of this L-series.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: L = E.lseries()
            sage: L._repr_()
            'Complex L-series of the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field'
        """
        return "Complex L-series of the %s"%self.__E

    def dokchitser(self, prec=53,
                   max_imaginary_part=0,
                   max_asymp_coeffs=40,
                   algorithm='gp'):
        r"""
        Return interface to Tim Dokchitser's program for computing
        with the L-series of this elliptic curve; this provides a way
        to compute Taylor expansions and higher derivatives of
        $L$-series.

        INPUT:
            prec -- integer (bits precision)
            max_imaginary_part -- real number
            max_asymp_coeffs -- integer
            algorithm -- string: 'gp' or 'magma'

        \note{If algorithm='magma', then the precision is in digits rather
        than bits and the object returned is a Magma L-series, which has
        different functionality from the SAGE L-series.}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: L = E.lseries().dokchitser()
            sage: L(2)
            0.381575408260711
            sage: L = E.lseries().dokchitser(algorithm='magma')         # optional - magma
            sage: L.Evaluate(2)                                         # optional - magma
            0.38157540826071121129371040958008663667709753398892116

        If the curve has too large a conductor, it isn't possible to
        compute with the L-series using this command.  Instead a
        RuntimeError is raised:
            sage: e = EllipticCurve([1,1,0,-63900,-1964465932632])
            sage: L = e.lseries().dokchitser(15)
            Traceback (most recent call last):
            ...
            RuntimeError: Unable to create L-series, due to precision or other limits in PARI.
        """
        if algorithm == 'magma':
            from sage.interfaces.all import magma
            return magma(self.__E).LSeries(Precision = prec)

        from sage.lfunctions.all import Dokchitser
        key = (prec, max_imaginary_part, max_asymp_coeffs)
        try:
            return self.__dokchitser[key]
        except KeyError:
            pass
        except AttributeError:
            self.__dokchitser = {}
        L = Dokchitser(conductor = self.__E.conductor(),
                       gammaV = [0,1],
                       weight = 2,
                       eps = self.__E.root_number(),
                       poles = [],
                       prec = prec)
        gp = L.gp()
        s = 'e = ellinit(%s);'%self.__E.minimal_model().a_invariants()
        s += 'a(k) = ellak(e, k);'
        L.init_coeffs('a(k)', 1, pari_precode = s,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        L.rename('Dokchitser L-function associated to %s'%self.__E)
        self.__dokchitser[key] = L
        return L

    def sympow(self, n, prec):
        r"""
        Return $L(\Sym^{(n)}(E, \text{edge}))$ to prec digits
        of precision.

        INPUT:
            n -- integer
            prec -- integer

        OUTPUT:
            string -- real number to prec digits of precision as a string.

        \note{Before using this function for the first time for
        a given $n$, you may have to type \code{sympow('-new_data <n>')},
        where \code{<n>} is replaced by your value of $n$.  This
        command takes a long time to run.}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: a = E.lseries().sympow(2,16)   # optional - requires precomputing "sympow('-new_data 2')"
            sage: a      # optional
            '2.492262044273650E+00'
            sage: RR(a)                      # optional
            2.49226204427365
        """
        from sage.lfunctions.sympow import sympow
        return sympow.L(self.__E, n, prec)

    def sympow_derivs(self, n, prec, d):
        r"""
        Return $0$th to $d$th derivatives of $L(\Sym^{(n)}(E,
        \text{edge}))$ to prec digits of precision.

        INPUT:
            n -- integer
            prec -- integer
            d -- integer

        OUTPUT:
            a string, exactly as output by sympow

        \note{To use this function you may have to run a few commands
        like \code{sympow('-new_data 1d2')}, each which takes a few
        minutes.  If this function fails it will indicate what
        commands have to be run.}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: print E.lseries().sympow_derivs(1,16,2)      # optional -- requires precomputing "sympow('-new_data 2')"
            sympow 1.018 RELEASE  (c) Mark Watkins --- see README and COPYING for details
            Minimal model of curve  is [0,0,1,-1,0]
            At 37: Inertia Group is  C1 MULTIPLICATIVE REDUCTION
            Conductor is 37
            sp 1: Conductor at 37 is 1+0, root number is 1
            sp 1: Euler factor at 37 is 1+1*x
            1st sym power conductor is 37, global root number is -1
            NT 1d0: 35
            NT 1d1: 32
            NT 1d2: 28
            Maximal number of terms is 35
            Done with small primes 1049
            Computed:  1d0  1d1  1d2
            Checked out:  1d1
             1n0: 3.837774351482055E-01
             1w0: 3.777214305638848E-01
             1n1: 3.059997738340522E-01
             1w1: 3.059997738340524E-01
             1n2: 1.519054910249753E-01
             1w2: 1.545605024269432E-01
        """
        from sage.lfunctions.sympow import sympow
        return sympow.Lderivs(self.__E, n, prec, d)

    def zeros(self, n):
        """
        Return the imaginary parts of the first $n$ nontrivial zeros
        on the critical line of the L-function in the upper half
        plane, as 32-bit reals.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.lseries().zeros(2)
            [0.000000000, 5.00317001]

            sage: a = E.lseries().zeros(20)             # long time
            sage: point([(1,x) for x in a]).save()    # graph  (long time)

        AUTHOR:
            -- Uses Rubinstein's L-functions calculator.
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.zeros(n, L=self.__E)

    def zeros_in_interval(self, x, y, stepsize):
        r"""
        Return the imaginary parts of (most of) the nontrivial zeros
        on the critical line $\Re(s)=1$ with positive imaginary part
        between $x$ and $y$, along with a technical quantity for each.

        INPUT:
            x, y, stepsize -- positive floating point numbers

        OUTPUT:
            list of pairs (zero, S(T)).

        Rubinstein writes: The first column outputs the imaginary part
        of the zero, the second column a quantity related to S(T) (it
        increases roughly by 2 whenever a sign change, i.e. pair of
        zeros, is missed). Higher up the critical strip you should use
        a smaller stepsize so as not to miss zeros.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.lseries().zeros_in_interval(6, 10, 0.1)      # long
            [(6.87039122, 0.248922780), (8.01433081, -0.140168533), (9.93309835, -0.129943029)]
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.zeros_in_interval(x, y, stepsize, L=self.__E)

    def values_along_line(self, s0, s1, number_samples):
        """
        Return values of $L(E, s)$ at \code{number_samples}
        equally-spaced sample points along the line from $s_0$ to
        $s_1$ in the complex plane.

        \note{The L-series is normalized so that the center of the
        critical strip is 1.}

        INPUT:
            s0, s1 -- complex numbers
            number_samples -- integer

        OUTPUT:
            list -- list of pairs (s, zeta(s)), where the s are
                    equally spaced sampled points on the line from
                    s0 to s1.

        EXAMPLES:
            sage: I = CC.0
            sage: E = EllipticCurve('37a')
            sage: E.lseries().values_along_line(1, 0.5+20*I, 5)     # long time and slightly random output
            [(0.500000000, 0), (0.400000000 + 4.00000000*I, 3.31920245 - 2.60028054*I), (0.300000000 + 8.00000000*I, -0.886341185 - 0.422640337*I), (0.200000000 + 12.0000000*I, -3.50558936 - 0.108531690*I), (0.100000000 + 16.0000000*I, -3.87043288 - 1.88049411*I)]
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.values_along_line(s0-RationalField()('1/2'),
                                       s1-RationalField()('1/2'),
                                       number_samples, L=self.__E)

    def twist_values(self, s, dmin, dmax):
        r"""
        Return values of $L(E, s, \chi_d)$ for each quadratic
        character $\chi_d$ for $d_{\min} \leq d \leq d_{\max}$.

        \note{The L-series is normalized so that the center of the
        critical strip is 1.}

        INPUT:
            s -- complex numbers
            dmin -- integer
            dmax -- integer

        OUTPUT:
            list -- list of pairs (d, L(E, s,chi_d))

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.lseries().twist_values(1, -12, -4)    # slightly random output depending on architecture
            [(-11, 1.4782434171), (-8, 0), (-7, 1.8530761916), (-4, 2.4513893817)]
            sage: F = E.quadratic_twist(-8)
            sage: F.rank()
            1
            sage: F = E.quadratic_twist(-7)
            sage: F.rank()
            0
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.twist_values(s - RationalField()('1/2'), dmin, dmax, L=self.__E)

    def twist_zeros(self, n, dmin, dmax):
        r"""
        Return first $n$ real parts of nontrivial zeros of
        $L(E,s,\chi_d)$ for each quadratic character $\chi_d$ with
        $d_{\min} \leq d \leq d_{\max}$.

        \note{The L-series is normalized so that the center of the
        critical strip is 1.}

        INPUT:
            n -- integer
            dmin -- integer
            dmax -- integer

        OUTPUT:
            dict -- keys are the discriminants $d$, and
                    values are list of corresponding zeros.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.lseries().twist_zeros(3, -4, -3)         # long
            {-4: [1.60813783, 2.96144840, 3.89751747], -3: [2.06170900, 3.48216881, 4.45853219]}
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.twist_zeros(n, dmin, dmax, L=self.__E)

    def at1(self, k=0):
        r"""
        Compute $L(E,1)$ using $k$ terms of the series for $L(E,1)$ as
        explained on page 406 of Henri Cohen's book"A Course in Computational
        Algebraic Number Theory".  If the argument $k$ is not specified,
        then it defaults to $\sqrt(N)$, where $N$ is the conductor.

        The real precision used in each step of the computation is the
        precision of machine floats.

        INPUT:
            k -- (optional) an integer, defaults to sqrt(N).

        OUTPUT:
            float -- L(E,1)
            float -- a bound on the error in the approximation; this
                     is a proveably correct upper bound on the sum
                     of the tail end of the series used to compute L(E,1).

        This function is disjoint from the PARI \code{elllseries}
        command, which is for a similar purpose.  To use that command
        (via the PARI C library), simply type
                \code{E.pari_mincurve().elllseries(1)}

        ALGORITHM:
        \begin{enumerate}
            \item Compute the root number eps.  If it is -1, return 0.

            \item Compute the Fourier coefficients a_n, for n up to and
               including k.

            \item Compute the sum
            $$
                 2 * sum_{n=1}^{k} (a_n / n) * exp(-2*pi*n/Sqrt(N)),
            $$
               where N is the conductor of E.

            \item Compute a bound on the tail end of the series, which is
            $$
                 2 * e^(-2 * pi * (k+1) / sqrt(N)) / (1 - e^(-2*pi/sqrt(N))).
            $$
               For a proof see [Grigov-Jorza-Patrascu-Patrikis-Stein].
        \end{enumerate}

        EXAMPLES:
            sage: E = EllipticCurve('37b')
            sage: E.lseries().at1(100)
            (0.725681061936153, 1.52437502288743e-45)
        """
        if self.__E.root_number() == -1:
            return 0
        sqrtN = float(self.__E.conductor().sqrt())
        k = int(k)
        if k == 0: k = int(ceil(sqrtN))
        an = self.__E.anlist(k)           # list of SAGE ints
        # Compute z = e^(-2pi/sqrt(N))
        pi = 3.14159265358979323846
        z = exp(-2*pi/sqrtN)
        zpow = z
        s = 0.0
        for n in xrange(1,k+1):
            s += (zpow * float(an[n]))/n
            zpow *= z

        error = 2*zpow / (1 - z)

        return R(2*s), R(error)

    def deriv_at1(self, k=0):
        r"""
        Compute $L'(E,1)$ using$ k$ terms of the series for $L'(E,1)$.

        The algorithm used is from page 406 of Henri Cohen's book ``A
        Course in Computational Algebraic Number Theory.''

        The real precision of the computation is the precision of
        Python floats.

        INPUT:
            k -- int; number of terms of the series

        OUTPUT:
            real number -- an approximation for L'(E,1)
            real number -- a bound on the error in the approximation

        ALGORITHM:
        \begin{enumerate}
            \item Compute the root number eps.  If it is 1, return 0.

            \item Compute the Fourier coefficients $a_n$, for $n$ up to and
                  including $k$.

            \item Compute the sum
               $$
                 2 * \sum_{n=1}^{k} (a_n / n) * E_1(2 \pi n/\sqrt{N}),
               $$
               where $N$ is the conductor of $E$, and $E_1$ is the
               exponential integral function.

            \item Compute a bound on the tail end of the series, which is
               $$
                 2 * e^{-2 \pi (k+1) / \sqrt{N}} / (1 - e^{-2 \ pi/\sqrt{N}}).
               $$
               For a proof see [Grigorov-Jorza-Patrascu-Patrikis-Stein].  This
               is exactly the same as the bound for the approximation to
               $L(E,1)$ produced by \code{E.lseries().at1}.
        \end{enumerate}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.lseries().deriv_at1()
            (0.305986660899342, 0.000800351433106958)
            sage: E.lseries().deriv_at1(100)
            (0.305999773834879, 1.52437502288740e-45)
            sage: E.lseries().deriv_at1(1000)
            (0.305999773834879, 0.000000000000000)
        """
        if self.__E.root_number() == 1: return 0
        k = int(k)
        sqrtN = float(self.__E.conductor().sqrt())
        if k == 0: k = int(ceil(sqrtN))
        an = self.__E.anlist(k)           # list of C ints
        # Compute z = e^(-2pi/sqrt(N))
        pi = 3.14159265358979323846
        v = transcendental.exponential_integral_1(2*pi/sqrtN, k)
        L = 2*float(sum([ (v[n-1] * an[n])/n for n in xrange(1,k+1)]))
        error = 2*exp(-2*pi*(k+1)/sqrtN)/(1-exp(-2*pi/sqrtN))
        return R(L), R(error)

    def __call__(self, s):
        r"""
        Returns the value of the L-series of the elliptic curve E at s, where s
        must be a real number.

        Use self.extended for s complex.

        \note{If the conductor of the curve is large, say $>10^{12}$,
        then this function will take a very long time, since it uses
        an $O(\sqrt{N})$ algorithm.}

        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: L = E.lseries()
            sage: L(1)
            0
            sage: L(1.1)
            0.285491007678148
            sage: L(1.1 + I)
            0.174851377216615 + 0.816965038124457*I
        """
        return self.dokchitser()(s)

    #def extended(self, s, prec):
    #    r"""
    #    Returns the value of the L-series of the elliptic curve E at s
    #    can be any complex number using prec terms of the power series
    #    expansion.
    #
    #
    #    WARNING: This may be slow.  Consider using \code{dokchitser()}
    #    instead.
    #
    #    INPUT:
    #        s -- complex number
    #        prec -- integer
    #
    #    EXAMPLES:
    #        sage: E = EllipticCurve('389a')
    #        sage: E.lseries().extended(1 + I, 50)
    #        -0.638409959099589 + 0.715495262192901*I
    #        sage: E.lseries().extended(1 + 0.1*I, 50)
    #        -0.00761216538818315 + 0.000434885704670107*I
    #
    #    NOTE: You might also want to use Tim Dokchitser's
    #    L-function calculator, which is available by typing
    #    L = E.lseries().dokchitser(), then evaluating L.  It
    #    gives the same information but is sometimes much faster.
    #
    #    """
    #    try:
    #        s = C(s)
    #    except TypeError:
    #        raise TypeError, "Input argument %s must be coercible to a complex number"%s
    #    prec = int(prec)
    #    if abs(s.imag()) < R(0.0000000000001):
    #        return self(s.real())
    #    N = self.__E.conductor()
    #    pi = R(constants.pi)
    #    Gamma = transcendental.gamma
    #    Gamma_inc = transcendental.gamma_inc
    #    a = self.__E.anlist(prec)
    #    eps = self.__E.root_number()
    #    sqrtN = float(N.sqrt())
    #    def F(n, t):
    #        return Gamma_inc(t+1, 2*pi*n/sqrtN) * C(sqrtN/(2*pi*n))**(t+1)
    #    return C(N)**(-s/2) * C(2*pi)**s * Gamma(s)**(-1)\
    #           * sum([a[n]*(F(n,s-1) + eps*F(n,1-s)) for n in xrange(1,prec+1)])

    def L1_vanishes(self):
        """
        Returns whether or not L(E,1) = 0. The result is provably
        correct if the Manin constant of the associated optimal
        quotient is <= 2.  This hypothesis on the Manin constant
        is true for all curves of conductor <= 40000 (by Cremona) and
        all semistable curves (i.e., squarefree conductor).

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.lseries().L1_vanishes()
            False
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.lseries().L1_vanishes()
            False
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.lseries().L1_vanishes()
            True
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.lseries().L1_vanishes()
            True
            sage: E = EllipticCurve([0, 0, 1, -38, 90])     # 361A (CM curve))
            sage: E.lseries().L1_vanishes()
            True
            sage: E = EllipticCurve([0,-1,1,-2,-1])         # 141C (13-isogeny)
            sage: E.lseries().L1_vanishes()
            False

        WARNING: It's conceivable that machine floats are not large
        enough precision for the computation; if this could be the
        case a RuntimeError is raised.  The curve's real period would
        have to be very small for this to occur.

        ALGORITHM: Compute the root number.  If it is -1 then L(E,s)
        vanishes to odd order at 1, hence vanishes.  If it is +1, use
        a result about modular symbols and Mazur's "Rational Isogenies"
        paper to determine a provably correct bound (assuming Manin
        constant is <= 2) so that we can determine whether L(E,1) = 0.

        AUTHOR: William Stein, 2005-04-20.
        """
        return self.L_ratio() == 0

    def L_ratio(self):
        r"""
        Returns the ratio $L(E,1)/\Omega$ as an exact rational
        number. The result is \emph{provably} correct if the Manin
        constant of the associated optimal quotient is $\leq 2$.  This
        hypothesis on the Manin constant is true for all semistable
        curves (i.e., squarefree conductor), by a theorem of Mazur
        from his \emph{Rational Isogenies of Prime Degree} paper.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.lseries().L_ratio()
            1/5
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.lseries().L_ratio()
            1/25
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.lseries().L_ratio()
            0
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.lseries().L_ratio()
            0
            sage: E = EllipticCurve([0, 0, 1, -38, 90])     # 361A (CM curve))
            sage: E.lseries().L_ratio()
            0
            sage: E = EllipticCurve([0,-1,1,-2,-1])         # 141C (13-isogeny)
            sage: E.lseries().L_ratio()
            1
            sage: E = EllipticCurve(RationalField(), [1, 0, 0, 1/24624, 1/886464])
            sage: E.lseries().L_ratio()
            2

            # See trac #3651:
            sage: EllipticCurve([0,0,0,-193^2,0]).sha().an()
            4

        WARNING: It's conceivable that machine floats are not large
        enough precision for the computation; if this could be the
        case a RuntimeError is raised.  The curve's real period would
        have to be very small for this to occur.

        ALGORITHM: Compute the root number.  If it is -1 then L(E,s)
        vanishes to odd order at 1, hence vanishes.  If it is +1, use
        a result about modular symbols and Mazur's "Rational Isogenies"
        paper to determine a provably correct bound (assuming Manin
        constant is <= 2) so that we can determine whether L(E,1) = 0.

        AUTHOR: William Stein, 2005-04-20.
        """
        try:
            return self.__lratio
        except AttributeError:
            pass

        if not self.__E.is_minimal():
            self.__lratio = self.__E.minimal_model().lseries().L_ratio()
            return self.__lratio

        if self.__E.root_number() == -1:
            self.__lratio = Q(0)
            return self.__lratio

        # Even root number.  Decide if L(E,1) = 0.  If E is a modular
        # *OPTIMAL* quotient of J_0(N) elliptic curve, we know that T *
        # L(E,1)/omega is an integer n, where T is the order of the
        # image of the rational torsion point (0)-(oo) in E(Q), and
        # omega is the least real Neron period.  (This is proved in my
        # Ph.D. thesis, but is probably well known.)  We can easily
        # compute omega to very high precision using AGM.  So to prove
        # that L(E,1) = 0 we compute T/omega * L(E,1) to sufficient
        # precision to determine it as an integer.  If eps is the
        # error in computation of L(E,1), then the error in computing
        # the product is (2T/Omega_E) * eps, and we need this to be
        # less than 0.5, i.e.,
        #          (2T/Omega_E) * eps < 0.5,
        # so
        #          eps < 0.5 * Omega_E / (2T) = Omega_E / (4*T).
        #
        # Since in general E need not be optimal, we have to choose
        # eps = Omega_E/(8*t*B), where t is the exponent of E(Q)_tor,
        # and is a multiple of the degree of an isogeny between E
        # and the optimal curve.
        #
        # NOTES: We *do* have to worry about the Manin constant, since
        # we are using the Neron model to compute omega, not the
        # newform.  My theorem replaces the omega above by omega/c,
        # where c is the Manin constant, and the bound must be
        # correspondingly smaller.  If the level is square free, then
        # the Manin constant is 1 or 2, so there's no problem (since
        # we took 8 instead of 4 in the denominator).  If the level
        # is divisible by a square, then the Manin constant could
        # be a divisible by an arbitrary power of that prime, except
        # that Edixhoven claims the primes that appear are <= 7.

        t = self.__E.torsion_subgroup().order()
        omega = self.__E.period_lattice().basis()[0]
        d = self.__E._multiple_of_degree_of_isogeny_to_optimal_curve()
        C = 8*d*t
        eps = omega / C
        # coercion of 10**(-15) to our real field is needed to
        # make unambiguous comparison
        if eps < R(10**(-15)):  # liberal bound on precision of float
            raise RuntimeError, "Insufficient machine precision (=%s) for computation."%eps
        sqrtN = 2*int(sqrt(self.__E.conductor()))
        k = sqrtN + 10
        while True:
            L1, error_bound = self.at1(k)
            if error_bound < eps:
                n = int(round(L1*C/omega))
                quo = Q(n) / Q(C)
                self.__lratio = quo / self.__E.real_components()
                return self.__lratio
            k += sqrtN
            misc.verbose("Increasing precision to %s terms."%k)
