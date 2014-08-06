r"""
Class file for computing sums over zeros of motivic L-functions.
All computations are done to double precision.

AUTHORS:

- Simon Spicer (2014-08): first version

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

from sage.structure.sage_object cimport SageObject
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
from sage.misc.all import verbose

cdef extern from "<math.h>":
    double c_exp "exp"(double)
    double c_log "log"(double)
    double c_cos "cos"(double)
    double c_acos "acos"(double)
    double c_sqrt "sqrt"(double)

cdef class LFunctionZeroSum_abstract(SageObject):
    """
    Abstract class for computing certain sums over zeros of a motivic L-function
    without having to determine the zeros themselves
    """
    cdef _N
    cdef _k
    cdef _C0
    cdef _C1
    cdef _pi
    cdef _euler_gamma

    def level(self):
        r"""
        Return the level of the form attached to self. If self was constructed
        from an elliptic curve, then this is equal to the conductor of 'E'.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: Z = LFunctionZeroSum(E)
            sage: Z.level()
            389

        """
        return self._N

    def weight(self):
        r"""
        Return the weight of the form attached to self. If self was constructed
        from an elliptic curve, then this is 2.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: Z = LFunctionZeroSum(E)
            sage: Z.weight()
            2

        """
        return self._k

    def C0(self, include_euler_gamma=True):
        r"""
        Return the constant term of the logarithmic derivative of the
        completed 'L'-function attached to self. This is equal to
            '-\eta + \log(N)/2 - \log(2\pi)'
        where '\eta' is the Euler-Mascheroni constant '= 0.5772\ldots'
        and 'N' is the level of the form attached to self.

        INPUT:

        - ``include_euler_gamma`` -- bool (default: True); if set to
          False, return the constant '\log(N)/2 - \log(2\pi)', i.e. do
          not subtract off the Euler-Mascheroni constant.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: Z = LFunctionZeroSum(E)
            sage: Z.C0()
            0.566696940498
            sage: Z.C0(include_euler_gamma=False)
            1.1439126054

        """
        # Computed at initialization
        if include_euler_gamma==False:
            return self._C1
        else:
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
        n = ZZ(2)
        while n <= num_terms:
            if n.is_prime_power():
                cn = self.cn(n)
                y += cn/F(n)**z
            n += 1

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

    def zerosum(self,Delta=1,tau=0,function="sincsquared_fast"):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*(\gamma-\tau),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x)' is an appropriate even continuous
        'L_2' function such that 'f(0)=1'.

        If '\tau=0', then as '\Delta' increases this sum limits from above to
        the analytic rank of the 'L'-function, as 'f(0) = 1' is counted with
        multiplicity 'r', and the other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter denoting the
          tightness of the zero sum.

        - ``tau`` -- real parameter (default: 0) denoting the offset of the sum
          to be computed. When 'tau=0' the sum will limit to the analytic rank
          of the 'L'-function as '\Delta' is increased. If 'tau' is the value
          of the imaginary part of a noncentral zero, the limit will be 1
          (assuming the zero is simple); otherwise the the limit will be 0.
          Currently only implemented for the sincsquared and cauchy functions;
          otherwise ignored.

        - ``function`` -- string (default: "sincsquared_fast") - the function
          'f(x)' as described above. Currently implemented options for 'f' are:

          - ``sincquared`` -- 'f(x) = \left(\frac{\sin(\pi*x)}{(\pi*x)}\right)^2'
          - ``gaussian``   -- 'f(x) = \exp(-x^2)'
          - ``sincquared_fast`` -- Same as "sincsquared", but implementation
            optimized for elliptic curve 'L'-functions with tau=0. self must
            be attached to an elliptic curve over QQ given by its global minimal
            model, otherwise the returned result will be incorrect.
          - ``cauchy`` -- f(x) = \frac{1}{1+x^2}; this is only computable to
            low precision, and only when Delta < 2.

        .. WARNING::

            Computation time is exponential in '\Delta', roughly doubling for
            every increase of 0.1 thereof. Using '\Delta=1' will yield a
            computation time of a few milliseconds; '\Delta=2' takes a few
            seconds, and '\Delta=3' takes upwards of an hour. Increase at your
            own risk beyond this!

        OUTPUT:

        A positive real number that bounds from above the number of zeros with
        imaginary part equal to tau. When tau=0 this is an upper bound for the
        'L'-function's analytic rank.

        .. SEEALSO::

            :meth:`~sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.analytic_rank_bound`
            for more documentation and examples on calling this method on elliptic curve
            'L'-functions.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); E.rank()
            2
            sage: Z = LFunctionZeroSum(E)
            sage: E.lseries().zeros(3)
            [0.000000000, 0.000000000, 2.87609907]
            sage: Z.zerosum(Delta=1,function="sincsquared_fast")
            2.0375000846
            sage: Z.zerosum(Delta=1,function="sincsquared")
            2.0375000846
            sage: Z.zerosum(Delta=1,tau=2.876,function='sincsquared')
            1.07555129565
            sage: Z.zerosum(Delta=1,tau=1.2,function='sincsquared')
            0.108315553775
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
            return self._zerosum_sincsquared(Delta=Delta,tau=tau)
        elif function=="gaussian":
            return self._zerosum_gaussian(Delta=Delta)
        elif function=="cauchy":
            return self._zerosum_cauchy(Delta=Delta,tau=tau)
        else:
            raise ValueError("Input function not recognized.")

    def _zerosum_sincsquared(self,Delta=1,tau=0):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*(\gamma-\tau)),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x) = \sin(\pi*x)/(\pi*x)'

        If '\tau=0', then as '\Delta' increases this sum limits from above to
        the analytic rank of the 'L'-function, as 'f(0) = 1' is counted with
        multiplicity 'r', and the other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter denoting the
          tightness of the zero sum.

        - ``tau`` -- real parameter (default: 0) denoting the offset of the sum
          to be computed. When 'tau=0' the sum will limit to the analytic rank
          of the 'L'-function as '\Delta' is increased. If 'tau' is the value
          of the imaginary part of a noncentral zero, the limit will be 1
          (assuming the zero is simple); otherwise the the limit will be 0.

        .. WARNING::

            Computation time is exponential in '\Delta', roughly doubling for
            every increase of 0.1 thereof. Using '\Delta=1' will yield a
            computation time of a few milliseconds; '\Delta=2' takes a few
            seconds, and '\Delta=3' takes upwards of an hour. Increase at your
            own risk beyond this!

        OUTPUT:

        A positive real number that bounds from above the number of zeros with
        imaginary part equal to tau. When tau=0 this is an upper bound for the
        'L'-function's analytic rank.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: E.lseries().zeros(2)
            [0.000000000, 5.00317001]

        E is a rank 1 curve; the lowest noncentral zero has imaginary part
        ~5.003. The zero sum with tau=0 indicates the probable existence of
        a zero at or very close to the the central point:

        ::

            sage: Z._zerosum_sincsquared(Delta=1,tau=0)
            1.01038406984

        The zero sum also detects a zero at or near 5.003, as expected:

        ::

            sage: Z._zerosum_sincsquared(Delta=1,tau=5.003)
            1.01681245469

        However, there is definitely no zero with imaginary part near 2.5,
        as the sum would have to be at least 1:

        ::

            sage: Z._zerosum_sincsquared(Delta=1,tau=2.5)
            0.0580582108065

        """

        npi = self._pi
        twopi = 2*npi
        eg = self._euler_gamma

        t = RDF(Delta*twopi)
        expt = RDF(exp(t))

        u = t*self.C0()

        # No offset: formulae are simpler
        if tau==0:
            w = RDF(npi**2/6-spence(1-RDF(1)/expt))

            y = RDF(0)
            n = int(1)
            while n < expt:
                cn  = self.cn(n)
                if cn!=0:
                    logn = log(RDF(n))
                    y += cn*(t-logn)
                    #print(n,cn*(t-logn))
                n += 1
        # When offset is nonzero, the digamma transform (w) must
        # be computed as an infinite sum
        else:
            tau = RDF(tau)
            cos_tau_t = (tau*t).cos()
            sin_tau_t = (tau*t).sin()
            w = RDF(0)
            for k in range(1,1001):
                a1 = tau**2/(k*(k**2+tau**2))
                a2 = (k**2-tau**2)/(k**2+tau**2)**2
                a3 = (2*k*tau)/(k**2+tau**2)**2

                w0 = a1*t + a2
                w0 -= (a2*cos_tau_t-a3*sin_tau_t)*exp(-k*t)
                w += w0

            y = RDF(0)
            n = int(1)
            while n < expt:
                cn  = self.cn(n)
                if cn!=0:
                    logn = log(RDF(n))
                    y += cn*(t-logn)*(tau*logn).cos()
                    #print(n,cn*(t-logn))
                n += 1

        #print(expt,t,u,w,y)
        #print

        return (u+w+y)*2/(t**2)

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

        t = RDF(Delta*npi*2)
        expt = RDF(exp(t))

        u = self.C0()

        w = RDF(0)
        for k in range(1,1001):
            w += RDF(1)/k-erfcx(Delta*k)*Deltasqrtpi

        y = RDF(0)
        n = int(1)

        # TO DO: Error analysis to make sure this bound is good enough to
        # avoid non-negligible trucation error
        while n < expt:
            cn  = self.cn(n)
            if cn != 0:
                logn = log(RDF(n))
                y += cn*exp(-(logn/(2*Delta))**2)
            n += 1

        # y is the truncation of an infinite sum, so we must add a value which
        # exceeds the max amount we could have left out.
        # WARNING: The value of 0.1 has been thumbsucked
        return RDF(u+w+y+0.1)/Deltasqrtpi

    def _zerosum_cauchy(self,Delta=1,tau=0,num_terms=None):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*(\gamma-\tau)),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x) = \frac{1}{1+x^2}'.

        If '\tau=0', then as '\Delta' increases this sum limits from above to
        the analytic rank of the 'L'-function, as 'f(0) = 1' is counted with
        multiplicity 'r', and the other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter denoting the
          tightness of the zero sum.

        - ``tau`` -- real parameter (default: 0) denoting the offset of the sum
          to be computed. When 'tau=0' the sum will limit to the analytic rank
          of the 'L'-function as '\Delta' is increased. If 'tau' is the value
          of the imaginary part of a noncentral zero, the limit will be 1
          (assuming the zero is simple); otherwise the the limit will be 0.

        - ``num_terms`` -- positive integer (default: None): the number of
          terms computed in the truncated Dirichlet series for the L-function
          attached to self. If left at None, this is set to
          ceil(exp(2*\pi*\Delta)), the same number of terms used in the other
          zero sum methods for this value of Delta.
          Increase num_terms to get more accuracy.

        .. WARNING::

            This value can only be provably computed when Delta < 2; an error
            will be thrown if a Delta value larger than 2 is supplied.
            Furthermore, beware that computation time is exponential in
            '\Delta', roughly doubling for every increase of 0.1 thereof.
            Using '\Delta=1' will yield a computation time of a few
            milliseconds, while '\Delta=2' takes a few seconds.

        OUTPUT:

        A positive real number that bounds from above the number of zeros with
        imaginary part equal to tau. When tau=0 this is an upper bound for the
        'L'-function's analytic rank.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: E.lseries().zeros(2)
            [6.36261389, 8.60353962]

        E is a rank zero curve; the lowest zero has imaginary part ~6.36. The
        zero sum with tau=0 indicates that there are no zeros at the central
        point (otherwise the returned value would be at least 1):

        ::

            sage: Z = LFunctionZeroSum(E)
            sage: Z._zerosum_cauchy(Delta=1,tau=0)
            0.970107398446

        The zero sum with tau=6.36 indicates there might be a zero in the
        vicinity:

        ::

            sage: Z._zerosum_cauchy(Delta=1,tau=6.36261389)
            2.18090462633

        However, there are no zeros with imaginary part close to 1.5:

        ::

            sage: Z._zerosum_cauchy(Delta=1,tau=1.5)
            0.982707203755

        Because of the weak convergence of the Dirichlet series close to the
        critical line, the bound will in general get *worse* for larger Delta.
        This can be mitigated somewhat by increasing the number of terms.

        ::

            sage: Z._zerosum_cauchy(Delta=1.5)
            12.9383525898
            sage: Z._zerosum_cauchy(Delta=1.5,num_terms=100000)
            10.3951839608

        An error will be thrown if a Delta value >= 2 is passed:

        ::
            sage: Z._zerosum_cauchy(Delta=2)
            Traceback (most recent call last):
            ...
            ValueError: Bound not provably computable for Delta >= 2
        """

        if Delta >= 2:
            raise ValueError("Bound not provably computable for Delta >= 2")
        Del = RDF(Delta)
        if num_terms is None:
            num_terms = int(exp(2*self._pi*Del))

        if tau==0:
            one = RDF(1)
            s = one/Del+one
            u,err = self.completed_logarithmic_derivative(s,num_terms)
        else:
            one = CDF(1)
            s = CDF(one/Del+one, tau)
            u,err = self.completed_logarithmic_derivative(s,num_terms)
            u = u.real()

        return (u+err)/Del

cdef class LFunctionZeroSum_EllipticCurve(LFunctionZeroSum_abstract):
    """
    Subclass for computing certain sums over zeros of an elliptic curve L-function
    without having to determine the zeros themselves.
    """
    cdef _E
    cdef _e

    def __init__(self,E,N=None):
        r"""
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
            self._N = E.conductor()
        # PARI minicurve for computing a_p coefficients
        self._e = E.pari_mincurve()

        self._pi = RDF(pi)
        self._euler_gamma = RDF(euler_gamma)

        # These constants feature in most (all?) sums over the L-function's zeros
        self._C1 = log(RDF(self._N))/2 - log(self._pi*2)
        self._C0 = self._C1 - self._euler_gamma

    def __repr__(self):
        r"""
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

    cdef double _sincsquared_summand_1(self,
                                       long long n,
                                       double t,
                                       int ap,
                                       double p,
                                       double logp,
                                       double thetap,
                                       double sqrtp,
                                       double logq,
                                       double thetaq,
                                       double sqrtq,
                                       double z):
        r"""
        Private cdef method to compute the logarithmic derivative
        summand for the sinc^2 sum at prime values for when
        n <= sqrt(bound), bound = exp(t)
        Called in self._zerosum_sincsquared_fast() method
        """
        ap = self._e.ellap(n)
        p = n
        sqrtp = c_sqrt(p)
        thetap = c_acos(ap/(2*sqrtp))
        logp = c_log(p)

        sqrtq = 1
        thetaq = 0
        logq = logp

        z = 0
        while logq < t:
            sqrtq *= sqrtp
            thetaq += thetap
            z += 2*c_cos(thetaq)*(t-logq)/sqrtq
            #print(c_exp(logq),-2*c_cos(thetaq)*(t-logq)/sqrtq*logp)
            logq += logp
        #print('a',n,t,-z*logp)
        return -z*logp

    cdef double _sincsquared_summand_2(self,
                                       long long n,
                                       double t,
                                       int ap,
                                       double p,
                                       double logp):
        r"""
        Private cdef method to compute the logarithmic derivative
        summand for the sinc^2 sum at prime values for when
        sqrt(bound) < n < bound, bound = exp(t)
        Called in self._zerosum_sincsquared_fast() method
        """
        ap = self._e.ellap(n)
        p = n
        logp = c_log(p)
        #print('b',n,t,-(t-logp)*(logp/p)*ap)
        return -(t-logp)*(logp/p)*ap

    cpdef _zerosum_sincsquared_fast(self,Delta=1,bad_primes=None):
        """
        A faster, more intelligent implementation of self._zerosum_sincsquared().

        .. NOTE::

            This will only produce correct output if self._E is given by its
            global minimal model, i.e. if self._E.is_minimal()==True.

        INPUT:

        - ``Delta`` -- positive real parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to self.
        - ``bad_primes`` -- (default: None) If not None, a list of primes dividing
          the level of the form attached to self. This is passable so that this
          method can be run on curves whose conductor is large enough to warrant
          precomputing bad primes.

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to self from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum_sincsquared`
            for the more general but slower version of this method.

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: Z = LFunctionZeroSum(E)
            sage: print(E.rank(),Z._zerosum_sincsquared_fast(Delta=1))
            (1, 1.01038406984)
            sage: E = EllipticCurve('121a')
            sage: Z = LFunctionZeroSum(E);
            sage: print(E.rank(),Z._zerosum_sincsquared_fast(Delta=1.5))
            (0, 0.0104712060087)
        """
        # If Delta>6.619, then we will most likely get overflow: some ap values
        # will be too large to fit into a c int
        if Delta > 6.619:
            raise ValueError("Delta value too large; will result in overflow")

        cdef double npi = self._pi
        cdef double twopi = npi*2
        cdef double eg = self._euler_gamma

        cdef double t,u,w,y,z,expt,bound1,logp,logq
        cdef double thetap,thetaq,sqrtp,sqrtq,p,q
        cdef int ap,aq

        cdef long long n
        cdef double N_double = self._N

        t = twopi*Delta
        expt = c_exp(t)

        u = t*(-eg + c_log(N_double)/2 - c_log(twopi))
        w = npi**2/6-spence(-expt**(-1)+1)

        y = 0
        # Do bad primes first. Add correct contributions and subtract
        # incorrect contribution, since we'll add them back later on.
        if bad_primes is None:
            bad_primes = self._N.prime_divisors()
        bad_primes = [prime for prime in bad_primes if prime<expt]
        #print(expt,bad_primes)
        for prime in bad_primes:
            n = prime
            ap = self._e.ellap(n)
            p = n
            sqrtp = c_sqrt(p)
            thetap = c_acos(ap/(2*sqrtp))
            logp = c_log(p)

            q = 1
            sqrtq = 1
            aq = 1
            thetaq = 0
            logq = logp

            z = 0
            while logq < t:
                q *= p
                sqrtq *= sqrtp
                aq *= ap
                thetaq += thetap
                # Actual value of this term
                z += (aq/q)*(t-logq)
                # Incorrect value of this term to be removed below
                z -= 2*c_cos(thetaq)*(t-logq)/sqrtq
                logq += logp
            y -= z*logp

        # Good prime case. Bad primes are treated as good primes, but their
        # contribution here is cancelled out above; this way we don't
        # have to check if each prime divides the level or not.

        # Must deal with n=2,3,5 separately
        for m in [2,3,5]:
            n = m
            if n<expt:
                y += self._sincsquared_summand_1(n,t,ap,p,logp,thetap,sqrtp,
                                                     logq,thetaq,sqrtq,z)
        # Now iteratonly only over those n that are 1 or 5 mod 6
        n = 11
        # First: those n that are <= sqrt(bound)
        bound1 = c_exp(t/2)
        while n <= bound1:
            if pari(n-4).isprime():
                y += self._sincsquared_summand_1(n-4,t,ap,p,logp,thetap,sqrtp,
                                                 logq,thetaq,sqrtq,z)
            if pari(n).isprime():
                y += self._sincsquared_summand_1(n,t,ap,p,logp,thetap,sqrtp,
                                                 logq,thetaq,sqrtq,z)
            n += 6
        # Unlucky split case where n-4 <= sqrt(bound) but n isn't
        if n-4 <= bound1 and n > bound1:
            if pari(n-4).isprime():
                y += self._sincsquared_summand_1(n-4,t,ap,p,logp,thetap,sqrtp,
                                                 logq,thetaq,sqrtq,z)
            if n <= expt and pari(n).isprime():
                y += self._sincsquared_summand_2(n,t,ap,p,logp)
            n += 6
        # Now sqrt(bound)< n < bound, so we don't need to consider higher
        # prime power logarithmic derivative coefficients
        while n <= expt:
            if pari(n-4).isprime():
                y += self._sincsquared_summand_2(n-4,t,ap,p,logp)
            if pari(n).isprime():
                y += self._sincsquared_summand_2(n,t,ap,p,logp)
            n += 6
        # Case where n-4 <= t but n isn't
        n = n-4
        if n <= expt and pari(n).isprime():
            y += self._sincsquared_summand_2(n,t,ap,p,logp)

        #print(expt,t,u,w,y)
        #print
        return RDF(2*(u+w+y)/(t**2))

    def analytic_rank_upper_bound(self,
                                  max_Delta=None,
                                  adaptive=True,
                                  root_number=True,
                                  bad_primes=None):
        """
        Return an upper bound for the analytic rank of the L-function
        'L_E(s)' attached to self, conditional on the Generalized Riemann
        Hypothesis, via computing the zero sum
            '\sum_{\gamma} f(\Delta*\gamma),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L(E,s)'
        along the critical strip, 'f(x) = (\sin(\pi x)/(\pi x))^2', and Delta
        is the tightness parameter whose maximum value is specified by
        max_Delta. This computation can be run on curves with very large
        conductor (so long as the conductor is known or quickly
        computable) when Delta is not too large (see below).
        Uses Bober's rank bounding method as described in [Bob13].

        INPUT:

        - ``max_Delta`` -- (default: None) If not None, a positive real value
          specifying the maximum Delta value used in the zero sum; larger
          values of Delta yield better bounds - but runtime is exponential in
          Delta. If left as None, Delta is set to
              '\min\{(\log(N+1000)/2-\log(2\pi)-\eta)/\pi\, 2.5}'
          where 'N' is the conductor of the curve attached to self, and '\eta'
          is the Euler-Mascheroni constant '= 0.5772\ldots'; the crossover
          point is at conductor ~8.3*10^8. For the former value, empirical
          results show that for about 99.7% of all curves the returned value
          is the actual analytic rank.

        - ``adaptive`` -- (default: True) Either True or False.
          - If True, the computation is first run with small and then
            successively larger Delta values up to max_Delta. If at any
            point the computed bound is 0 (or 1 when when root_number is -1
            or True), the computation halts and that value is returned;
            otherwise the minimum of the computed bounds is returned.
          - If False, the computation is run a single time with
            Delta=max_Delta, and the resulting bound returned.

        - ``root_number`` -- (default: True) One of True, False, 1 or -1:
          - If True, the root number of self is computed and used to (possibly)
            lower ther analytic rank estimate by 1.
          - If False, the above step is omitted.
          - If 1 or -1, this value is assumed to be the root number of self.
            This is passable so that rank estimation can be done for curves
            whose root number has been precomputed.

        - ``bad_primes`` -- (default: None) If not None, a list of the primes
          of bad reduction for the curve attached to self. This is passable
          so that rank estimation can be done for curves of large conductor
          whose bad primes have been precomputed.

        .. NOTE::

            Output will be incorrect if the incorrect root number is specified.

        .. WARNING::

            Zero sum computation time is exponential in the tightness parameter
            '\Delta', roughly doubling for every increase of 0.1 thereof.
            Using '\Delta=1' (and adaptive=False) will yield a runtime of a few
            milliseconds; '\Delta=2' takes a few seconds, and '\Delta=3' may
            take upwards of an hour. Increase beyond this at your own risk!

        OUTPUT:

        A non-negative integer greater than or equal to the analytic rank of
        self.

        .. NOTE::

            If you use set_verbose(1), extra information about the computation
            will be printed.

        .. SEEALSO::

            :func:`LFunctionZeroSum`
            :meth:`EllipticCurve.root_number`
            :func:`set_verbose`

        EXAMPLES:

        For most elliptic curves with small conductor the central zero(s)
        of 'L_E(s)' are fairly isolated, so small values of '\Delta'
        will yield tight rank estimates:

        ::

            sage: E = EllipticCurve('11a')
            sage: E.rank()
            0
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound(max_Delta=1)
            0

            sage: E = EllipticCurve([-39,123])
            sage: E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound(max_Delta=1)
            1

        This is especially true for elliptic curves with large rank:

        ::

            sage: E = elliptic_curves.rank(8)[0]; E
            Elliptic Curve defined by y^2 + y = x^3 - 23737*x + 960366 over Rational Field
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound(max_Delta=1)
            8

        However, some curves have 'L'-functions with low-lying zeroes, and for these
        larger values of '\Delta' must be used to get tight estimates:

        ::

            sage: E = EllipticCurve('974b1')
            sage: r = E.rank(); r
            0
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound(max_Delta=1,root_number=False)
            1
            sage: Z.analytic_rank_upper_bound(max_Delta=1.3,root_number=False)
            0

        Knowing the root number of E allows us to use smaller Delta values
        to get tight bounds, thus speeding up runtime considerably:

        ::

            sage: Z.analytic_rank_upper_bound(max_Delta=0.6,root_number=True)
            0

        REFERENCES:

        [Bob13] J.W. Bober. Conditionally bounding analytic ranks of elliptic curves.
        ANTS 10. http://msp.org/obs/2013/1-1/obs-v1-n1-p07-s.pdf

        """
        #Helper function: compute zero sum and apply parity if not False
        def run_computation(Delta,parity,bad_primes):
            verbose("Computing zero sum with Delta = %s"%Delta)
            bound = self._zerosum_sincsquared_fast(Delta=Delta,
                                                   bad_primes=bad_primes)
            verbose("Sum value is %s"%bound)
            bound = bound.floor()
            # parity is set to -1 when we're not taking root number into
            # account
            if parity==-1:
                verbose("Without invoking parity, rank bound is %s"%bound)
                return bound
            # parity is 0 if E has even analytic rank, and 1 if odd
            # analytic rank. The returned value must have the same parity
            # as the parity parameter.
            if bound%2!=parity:
                bound -= 1
            verbose("Invoking parity, rank bound is %s"%bound)
            return bound

        # Get/compute parity
        if root_number==False:
            verbose("Curve parity ignored.")
            parity = -1
        elif root_number==True:
            verbose("Computing curve parity...")
            parity = (1-self._e.ellrootno())//2
            verbose("Curve has parity %s."%parity)
        else:
            parity = (1-root_number)//2
            verbose("Parity set to %s."%parity)
        if parity==1:
            halt_bound = 1
            verbose("Computation will halt if at any point bound is <= 1.")
        else:
            halt_bound = 0
            verbose("Computation will halt if at any point bound is 0.")

        # Compute max_Delta if necessary
        if max_Delta is None:
            verbose("Computing maximum Delta value")
            pi, eg = self._pi, self._euler_gamma
            #1000 is arbitrary - increases Delta for small N
            max_Delta = (log(RDF(self._N+1000))/2-log(2*pi)-eg)/pi
            if max_Delta > 2.5:
                max_Delta = 2.5
                verbose("Computed max Delta value too big; setting to 2.5")
            else:
                verbose("Maximum Delta value to be used set at %s"%max_Delta)
        else:
            verbose("Maximum Delta value to be used set at %s"%max_Delta)

        # When max_Delta <= 1 it's not worth running the computation
        # multiple times, as it's so quick anyway
        if adaptive==False or max_Delta<=1:
            return run_computation(max_Delta,parity,bad_primes)
        else:
            bound_list = []
            # Find starting value. This loop won't ever take long,
            # since max_Delta is never > 7.
            Delta = max_Delta
            while Delta > 1:
                Delta -= 0.2
            # Now go back up the sequence of Deltas
            while Delta <= max_Delta:
                bound = run_computation(Delta,parity,bad_primes)
                if bound <= halt_bound:
                    verbose("computed bound <= halt_bound, so halting")
                    return bound
                else:
                    bound_list.append(bound)
                # Incrementing Delta by 0.2 each step means runtime
                # will increase by a factor of ~3.7
                Delta += 0.2

            # Since the zero sum is not strictly decreasing in Delta,
            # the last value is not necessarily the smallest
            smallest_bound = min(bound_list)
            verbose("Smallest bound computed is %s"%smallest_bound)
            return smallest_bound

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
