"""
Class file for computing sums over zeros of motivic L-functions.
All computations are done to double precision.

AUTHORS:

- Simon Spicer (2014-06): first version

"""

##############################################################################
#       Copyright (C) 2014 Simon Spicer <mlungu@uw.edu>
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
##############################################################################

from sage.structure.sage_object import SageObject
from scipy.special import erfcx, spence, psi
from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.infinity import PlusInfinity
from sage.rings.arith import prime_powers
from sage.functions.log import log, exp
from sage.functions.other import real, imag
from sage.symbolic.constants import pi, euler_gamma
from sage.libs.pari.all import pari

class LFunctionZeroSum_abstract(SageObject):
    """
    Abstract class for computing certain sums over zeros of a motivic L-function
    without having to determine the zeros themselves
    """

    def level(self):
        """
        Return the level of the form attached to self. If self was constructed
        from an elliptic curve, then this is equal to the conductor of E.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: Z = LFunctionZeroSum(E)
            sage: Z.level()
            389

        """
        return self._N

    def weight(self):
        """
        Return the weight of the form attached to self. If self was constructed
        from an elliptic curve, then this is 2.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: Z = LFunctionZeroSum(E)
            sage: Z.weight()
            2

        """
        return self._k

    def C0(self):
        r"""
        Return the constant term of the logarithmic derivative of the
        completed 'L'-function attached to self. This is equal to
            '-\eta + \log(N)/2 - \log(2\pi)'
        where '\eta' is the Euler-Mascheroni constant '= 0.5772\ldots'
        and 'N' is the level of the form attached to self.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: Z = LFunctionZeroSum(E)
            sage: Z.C0()
            0.566696940498

        """
        # Defined at initialization
        return self._C0

    def cnlist(self,n,python_floats=False):
        r"""
        Return a list of Dirichlet coefficient of the logarithmic
        derivative of the L-function attached to self, shifted so that
        the critical line lies on the imaginary axis, up to and
        including n. The i-th element of the return list is a[i].

        INPUT:

        - ``n`` -- non-negative integer

        - ``python_floats`` -- bool (default: False); if True return a list of
          Python floats instead of Sage Real Double Field elements.

        OUTPUT:

        A list of real numbers

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_EllipticCurve.cn`

        .. TODO::

            Speed this up; make more efficient

        EXAMPLES:

            sage: E = EllipticCurve('11a')
            sage: Z = LFunctionZeroSum(E)
            sage: cnlist = Z.cnlist(11)
            sage: for n in range(12): print(n,cnlist[n])
            (0, 0.0)
            (1, 0.0)
            (2, 0.69314718056)
            (3, 0.366204096223)
            (4, 0.0)
            (5, -0.321887582487)
            (6, 0.0)
            (7, 0.555974328302)
            (8, -0.34657359028)
            (9, 0.610340160371)
            (10, 0.0)
            (11, -0.217990479345)
        """
        if not python_floats:
            return [self.cn(i) for i in xrange(n+1)]
        else:
            return [float(self.cn(i)) for i in xrange(n+1)]

    def digamma(self,s,include_constant_term=True):
        r"""
        Return the digamma function on the complex input 's', given by
            '\digamma(s) = -\eta + \sum_{k=1)^{\infty} frac{s-1}{k(k+s-1)}'
        where '\eta' is the Euler-Mascheroni constant '=0.5772156649\ldots'.
        This function is needed in the computing the logarithmic derivative
        of the 'L'-function attached to self.

        INPUT:

        - ``s`` -- A complex number
        - ``include_constant_term`` -- (default: True) boolean; if set False,
          only the value of the sum over 'k' is returned without subtracting
          off the Euler-Mascheroni constant, i.e. the returned value is
          equal to '\sum_{k=1)^{\infty} frac{s-1}{k(k+s-1)}'.

        OUTPUT:

        A real double precision number if the input is real and not a negative
        integer; Infinity if the input is a negative integer, and a complex
        number otherwise.

        EXAMPLES::

            sage: Z = LFunctionZeroSum(EllipticCurve('37a'))
            sage: Z.digamma(3.2)
            0.998838891287
            sage: Z.digamma(3.2,include_constant_term=False)
            1.57605455619
            sage: Z.digamma(1+I)
            0.0946503206225 + 1.07667404747*I
            sage: Z.digamma(-2)
            +Infinity

        Evaluating the sum without the constant term at the positive integers n
        returns the (n-1)th harmonic number:

        ::

            sage: Z.digamma(3,include_constant_term=False)
            1.5
            sage: Z.digamma(6,include_constant_term=False)
            2.28333333333
        """
        if real(s)<0 and imag(s)==0:
            try:
                z = ZZ(s)
                return PlusInfinity()
            except:
                pass

        if imag(s)==0:
            F = RDF
        else:
            F = CDF
        # Cheating: SciPy already has this function implemented
        z = F(psi(F(s)))
        if include_constant_term:
            return z
        else:
            return z + self._euler_gamma

    def logarithmic_derivative(self,s,num_terms=10000,as_interval=False):
        r"""
        Compute the value of the logarithmic derivative '\frac{L^{\prime}}{L}'
        at the point s to *low* precision, where 'L' is the L-function
        attached to self.

        .. WARNING::

            The value is computed naively by evaluating the Dirichlet series
            for '\frac{L^{\prime}}{L}'; convergence is controlled by the
            distance of s from the critical strip '0.5<=\Re(s)<=1.5'.
            You may use this method to attempt to compute values inside the
            critical strip; however, results are then *not* guaranteed
            to be correct to any number of digits.

        INPUT:

        - ``s`` -- Real or complex value
        - ``num_terms`` -- (default: 10000) the maximum number of terms
          summed in the Dirichlet series.

        OUTPUT:

        A tuple (z,err), where z is the computed value, and err is an
        upper bound on the truncation error in this value introduced
        by truncating the Dirichlet sum.

        .. NOTE::

        For the default term cap of 10000, a value accurate to all 53
        bits of a double precision floating point number is only
        guaranteed when '|Re(s-1)|>4.58', although in practice inputs
        closer to the critical strip will still yield computed values
        close to the true value.

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.logarithmic_derivative(10)
            (5.64806674263e-05, 1.09741028598e-34)
            sage: Z.logarithmic_derivative(2.2)
            (0.575125706359, 0.024087912697)

        Increasing the number of terms should see the truncation error
        decrease:

        ::

            sage: Z.logarithmic_derivative(2.2,num_terms=50000) # long time
            (0.575157964506, 0.00898877551916)

        Attempting to compute values inside the critical strip
        gives infinite error:

        ::

            sage: Z.logarithmic_derivative(1.3)
            (5.44299441392, +Infinity)

        Complex inputs and inputs to the left of the critical strip
        are allowed:

        ::

            sage: Z.logarithmic_derivative(complex(3,-1))
            (0.0476454857805 + 0.1651383281*I, 6.5846713591e-06)
            sage: Z.logarithmic_derivative(complex(-3,-1.1))
            (-13.9084521732 + 2.59144309907*I, 2.71315847363e-14)

        The logarithmic derivative has poles at the negative integers:

        ::

            sage: Z.logarithmic_derivative(-3)
            (-Infinity, 2.71315847363e-14)
        """
        if imag(s)==0:
            F = RDF
        else:
            F = CDF
        # Inputs left of the critical line are handled via the functional
        # equation of the logarithmic derivative
        if real(s-1)<0:
            a = -2*self._C1-self.digamma(s)-self.digamma(2-s)
            b,err = self.logarithmic_derivative(2-s,num_terms=num_terms)
            return (a+b,err)

        z = s-1
        sigma = RDF(real(z))
        log2 = log(RDF(2))
        # Compute maximum possible Dirichlet series truncation error
        # When s is in the critical strip: no guaranteed precision
        if abs(sigma)<=0.5:
            err = PlusInfinity()
        else:
            a = RDF(sigma)-RDF(0.5)
            b = log(RDF(num_terms))*a
            err = (b+1)*exp(-b)/a**2

        y = F(0)
        for n in prime_powers(2,num_terms+1):
            cn = self.cn(n)
            y += cn/F(n)**z

        return (y,err)

    def completed_logarithmic_derivative(self,s,num_terms=10000):
        r"""
        Compute the value of the completed logarithmic derivative
            '\frac{\Lambda^{\prime}}{\Lambda}'
        at the point s to *low* precision, where
            '\Lambda = N^{s/2}(2\pi)^s \Gamma(s) L(s)'
        and 'L is the' L-function attached to self.

        .. WARNING::

            This is computed naively by evaluating the Dirichlet series
            for '\frac{L^{\prime}}{L}'; the convergence thereof is
            controlled by the distance of s from the critical strip
            '0.5<=\Re(s)<=1.5'.
            You may use this method to attempt to compute values inside the
            critical strip; however, results are then *not* guaranteed
            to be correct to any number of digits.

        INPUT:

        - ``s`` -- Real or complex value
        - ``num_terms`` -- (default: 10000) the maximum number of terms
          summed in the Dirichlet series.

        OUTPUT:

        A tuple (z,err), where z is the computed value, and err is an
        upper bound on the truncation error in this value introduced
        by truncating the Dirichlet sum.

        .. NOTE::

        For the default term cap of 10000, a value accurate to all 53
        bits of a double precision floating point number is only
        guaranteed when '|Re(s-1)|>4.58', although in practice inputs
        closer to the critical strip will still yield computed values
        close to the true value.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_EllipticCurve.logarithmic_derivative`

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.completed_logarithmic_derivative(3)
            (6.64372066048, 6.5846713591e-06)

        Complex values are handled. The function is odd about s=1, so
        the value at 2-s should be minus the value at s:

            sage: Z.completed_logarithmic_derivative(complex(-2.2,1))
            (-6.89808063313 + 0.225570153942*I, 5.62385304981e-11)
            sage: Z.completed_logarithmic_derivative(complex(4.2,-1))
            (6.89808063313 - 0.225570153942*I, 5.62385304981e-11)
        """
        if imag(s)==0:
            F = RDF
        else:
            F = CDF

        if real(s-1)>=0:
            Ls = self.logarithmic_derivative(s,num_terms)
            return (self._C1 + self.digamma(s) + Ls[0], Ls[1])
        else:
            Ls = self.logarithmic_derivative(2-s,num_terms)
            return (-self._C1 - self.digamma(2-s) - Ls[0], Ls[1])

    def zerosum(self,Delta=1,function="sincsquared_fast"):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*\gamma),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x)' is an appropriate continuous
        'L_2' function such that 'f(0)=1'.

        As '\Delta' increases this sum limits from above to the analytic rank
        of the form, as 'f(0) = 1' is counted with multiplicity 'r', and the
        other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to self.

        - ``function`` -- string (default: "sincsquared_fast") - the function
          'f(x)' as described above. Currently implemented options for 'f' are:

          - ``sincquared`` -- 'f(x) = \left(\frac{\sin(\pi*x)}{(\pi*x)}\right)^2'
          - ``gaussian``   -- 'f(x) = \exp(-x^2)'
          - ``sincquared_fast`` -- Same as "sincsquared", but faster
            implementation; however self must be attached to an elliptic curve
            over QQ given by its global minimal model, otherwise the returned
            result will be incorrect.
          - ``cauchy`` -- f(x) = \frac{1}{1+x^2}; this is only computable to
            low precion, and only when Delta < 2.

        .. WARNING::

            Computation time is exponential in '\Delta', roughly doubling for
            every increase of 0.1 thereof. Using '\Delta=1' will yield a
            computation time of a few milliseconds; '\Delta=2' takes a few
            seconds, and '\Delta=3' takes upwards of an hour. Increase at your
            own risk beyond this!

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to self from above.

        .. SEEALSO::

            :meth:`~sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.analytic_rank_bound`
            for more documentation and examples on calling this method on elliptic curve
            'L'-functions.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); E.rank()
            2
            sage: Z = LFunctionZeroSum(E)
            sage: Z.zerosum(Delta=1,function="sincsquared_fast")
            2.0375000846
            sage: Z.zerosum(Delta=1,function="sincsquared")
            2.0375000846
            sage: Z.zerosum(Delta=1,function="gaussian")
            2.05689042503

        """

        # If Delta>6.95, then exp(2*pi*Delta)>sys.maxint, so we get overflow
        # when summing over the logarithmic derivative coefficients
        if Delta > 6.95:
            raise ValueError("Delta value too large; will result in overflow")

        if function=="sincsquared_fast":
            return self._zerosum_sincsquared_fast(Delta=Delta)
        elif function=="sincsquared":
            return self._zerosum_sincsquared(Delta=Delta)
        elif function=="gaussian":
            return self._zerosum_gaussian(Delta=Delta)
        elif function=="cauchy":
            return self._zerosum_cauchy(Delta=Delta)
        else:
            raise ValueError("Input function not recognized.")

    def _zerosum_sincsquared(self,Delta=1):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*\gamma),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x) = \sin(\pi*x)/(\pi*x)'

        As '\Delta' increases this sum limits from above to the analytic rank
        of the form, as 'f(0) = 1' is counted with multiplicity 'r', and the
        other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to self.

        .. WARNING::

            Computation time is exponential in '\Delta', roughly doubling for
            every increase of 0.1 thereof. Using '\Delta=1' will yield a
            computation time of a few milliseconds; '\Delta=2' takes a few
            seconds, and '\Delta=3' takes upwards of an hour. Increase at your
            own risk beyond this!

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to self from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z._zerosum_sincsquared(Delta=1)
            1.01038406984

        """

        npi = self._pi
        twopi = 2*npi
        eg = self._euler_gamma

        t = RDF(Delta*twopi)
        expt = RDF(exp(t))

        u = t*(-eg + log(RDF(self._N))/2 - log(twopi))
        w = RDF(npi**2/6-spence(1-RDF(1)/expt))

        y = RDF(0)
        n = int(0)
        while n < expt:
            n += 1
            cn  = self.cn(n)
            if cn!=0:
                logn = log(RDF(n))
                y += cn*(t-logn)

        return 2*(u+w+y)/(t**2)

    def _zerosum_gaussian(self,Delta=1):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*\gamma),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x) = \exp(-x^2)'

        As '\Delta' increases this sum limits from above to the analytic rank
        of the form, as 'f(0) = 1' is counted with multiplicity 'r', and the
        other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to self.

        .. WARNING::

            Computation time is exponential in '\Delta', roughly doubling for
            every increase of 0.1 thereof. Using '\Delta=1' will yield a
            computation time of a few milliseconds; '\Delta=2' takes a few
            seconds, and '\Delta=3' takes upwards of an hour. Increase at your
            own risk beyond this!

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to self from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z._zerosum_gaussian(Delta=1)
            1.05639507734

        """

        npi = self._pi
        eg = self._euler_gamma
        Deltasqrtpi = Delta*npi.sqrt()

        u = -eg + log(RDF(self._N))/2 - log(npi*2)

        w = RDF(0)
        for k in range(1,1001):
            w += RDF(1)/k-erfcx(Delta*k)*Deltasqrtpi

        y = RDF(0)
        n = int(0)
        exp2piDelta = exp(2*npi*Delta)

        # TO DO: Error analysis to make sure this bound is good enough to
        # avoid non-negligible trucation error
        while n < exp2piDelta:
            n += 1
            cn  = self.cn(n)
            if cn != 0:
                logn = log(RDF(n))
                y += cn*exp(-(logn/(2*Delta))**2)

        # y is the truncation of an infinite sum, so we must add a value which
        # exceeds the max amount we could have left out.
        return RDF(u+w+y+0.1)/Deltasqrtpi

    def _zerosum_cauchy(self,Delta=1):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*\gamma),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x) = \frac{1}{1+x^2}'.

        As '\Delta' increases this sum limits from above to the analytic rank
        of the form, as 'f(0) = 1' is counted with multiplicity 'r', and the
        other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to self.

        .. WARNING::

            This value can only be computed when Delta < 2. An error will be
            thrown if a Delta value larger than 2 is supplied. Furthermore,
            beware that computation time is exponential in '\Delta', roughly
            doubling for every increase of 0.1 thereof. Using '\Delta=1' will
            yield a computation time of a few milliseconds, while '\Delta=2'
            takes a few seconds.

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to self from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

        """
        raise NotImplemedError("Method not yet implemented.")

class LFunctionZeroSum_EllipticCurve(LFunctionZeroSum_abstract):
    """
    Subclass for computing certain sums over zeros of an elliptic curve L-function
    without having to determine the zeros themselves
    """
    def __init__(self,E,N=None):
        """
        Initializes self.

        INPUT:

        - ``E`` -- An elliptic curve defined over the rational numbers

        - ``N`` -- (default: None) If not None, a positive integer equal to
          the conductor of E. This is passable so that rank estimation
          can be done for curves whose (large) conductor has been precomputed.

        EXAMPLES:

            sage: from sage.lfunctions.zero_sums import LFunctionZeroSum_EllipticCurve
            sage: E = EllipticCurve([1,0,0,3,-4])
            sage: Z = LFunctionZeroSum_EllipticCurve(E); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + x*y = x^3 + 3*x - 4 over Rational Field
            sage: E = EllipticCurve('5077a')
            sage: Z = LFunctionZeroSum_EllipticCurve(E); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
        """

        self._k = ZZ(2)
        self._E = E
        if N is not None:
            self._N = N
        else:
            self._N = self._E.conductor()
        # PARI minicurve for computing a_p coefficients
        self._e = E.pari_mincurve()

        self._pi = RDF(pi)
        self._euler_gamma = RDF(euler_gamma)

        # These constants feature in most (all?) sums over the L-function's zeros
        self._C1 = log(RDF(self._N))/2 - log(self._pi*2)
        self._C0 = self._C1 - self._euler_gamma

    def __repr__(self):
        """
        Representation of self.

        EXAMPLES::

            sage: Z = LFunctionZeroSum(EllipticCurve('37a')); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        s = "Zero sum estimator for L-function attached to "
        return s+str(self._E)

    def elliptic_curve(self):
        """
        Return the elliptic curve associated with self.

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.elliptic_curve()
            Elliptic Curve defined by y^2 = x^3 + 23*x + 100 over Rational Field

        """
        return self._E

    def lseries(self):
        """
        Return the 'L'-series associated with self.

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.lseries()
            Complex L-series of the Elliptic Curve defined by y^2 = x^3 + 23*x + 100 over Rational Field

        """
        return self._E.lseries()

    def cn(self,n):
        r"""
        Return the 'n'th Dirichlet coefficient of the logarithmic
        derivative of the L-function attached to self, shifted so that
        the critical line lies on the imaginary axis. This is zero if
        'n' is not a perfect prime power; and when 'n=p^e' it is
        '-a_p^e*log(p)/p^e', where 'a_p = p+1-\#{E(FF_p)}' is the
        trace of Frobenius at 'p'.

        INPUT:

        - ``n`` -- non-negative integer

        OUTPUT:

        A real number at most '\frac{log(n)}{\sqrt{n}}' in magnitude

        EXAMPLES::

        sage: E = EllipticCurve('11a')
        sage: Z = LFunctionZeroSum(E)
        sage: for n in range(12): print(n,Z.cn(n))
        (0, 0.0)
        (1, 0.0)
        (2, 0.69314718056)
        (3, 0.366204096223)
        (4, 0.0)
        (5, -0.321887582487)
        (6, 0.0)
        (7, 0.555974328302)
        (8, -0.34657359028)
        (9, 0.610340160371)
        (10, 0.0)
        (11, -0.217990479345)

        """
        n = ZZ(n)
        if n==0 or n==1:
            return RDF(0)
        if not n.is_prime_power():
            return RDF(0)

        n_float = RDF(n)
        if n.is_prime():
            logn = log(n_float)
            ap = self._E.ap(n)
            return -ap*logn/n_float
        else:
            p,e = n.perfect_power()
            ap = self._E.ap(p)
            logp = log(RDF(p))
            if p.divides(self._N):
                return - ap**e*logp/n_float
            c = CDF(ap,(4*p-ap**2).sqrt())/2
            aq = (2*(c**e).real()).round()
            return -aq*logp/n_float

    def _zerosum_sincsquared_fast(self,Delta=1):
        """
        A faster, more intelligent implementation of self._zerosum_sincsquared().

        .. NOTE::

            This will only produce correct output if self._E is given by its
            global minimal model, i.e. if self._E.is_minimal()==True.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to self.

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to self from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum_sincsquared`
            for the more general but slower version of this method.

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z._zerosum_sincsquared_fast(Delta=1)
            1.01038406984
        """

        npi = self._pi
        twopi = 2*npi
        eg = self._euler_gamma

        t = RDF(Delta*twopi)
        expt = exp(t)

        u = t*(-eg + log(RDF(self._N))/2 - log(twopi))
        w = RDF(npi**2/RDF(6)-spence(RDF(1)-RDF(1)/expt))

        y = RDF(0)
        bound1 = int(exp(t/2))
        bound2 = int(expt)
        n = int(2)
        while n <= bound1:
            if pari(n).isprime():
                logp = log(RDF(n))
                ap = RDF(self._e.ellap(n))
                p = RDF(n)

                z = (ap/p)*(t-logp)

                # The p^n coefficients are calculated differently
                # depending on whether p divides the level or not
                if self._N%n==0:
                    aq = ap**2
                    q = p**2
                    logq = logp*2
                    while logq < t:
                        z += (aq/q)*(t-logq)
                        logq += logp
                        q = q*p
                        aq = aq*ap
                else:
                    alpha_p = CDF(ap,(4*p-ap**2).sqrt())/2
                    alpha = alpha_p**2
                    q = p**2
                    logq = logp*2
                    while logq < t:
                        aq = RDF(round(2*alpha.real()))
                        z += (aq/q)*(t-logq)
                        logq += logp
                        q = q*p
                        alpha = alpha*alpha_p
                y -= z*logp
            n += 1
        while n <= bound2:
            if pari(n).isprime():
                logp = log(RDF(n))
                ap = RDF(self._e.ellap(n))
                y -= ap*logp/RDF(n)*(t-logp)
            n += 1

        return 2*(u+w+y)/(t**2)


def LFunctionZeroSum(X,*args,**kwds):
    """
    Constructor for the LFunctionZeroSum class.

    INPUT:

    - ''X'' - A motivic object. Currently only implemented for X = an elliptic curve
      over the rational numbers.

    OUTPUT:

    An LFunctionZeroSum object.

    EXAMPLES::

        sage: E = EllipticCurve('389a')
        sage: Z = LFunctionZeroSum(E); Z
        Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    TESTS::

        sage: E = EllipticCurve([1.2,3.8])
        sage: LFunctionZeroSum(E)
        Traceback (most recent call last):
        ...
        NotImplementedError: Currently only implemented for elliptic curves over QQ

        sage: f = Newforms(46)[0]
        sage: LFunctionZeroSum(f)
        Traceback (most recent call last):
        ...
        NotImplementedError: Currently only implemented for elliptic curves over QQ

    """

    # Here to avoid import recursion
    from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field

    if isinstance(X,EllipticCurve_rational_field):
        return LFunctionZeroSum_EllipticCurve(X,*args,**kwds)

    raise NotImplementedError("Currently only implemented for elliptic curves over QQ")
