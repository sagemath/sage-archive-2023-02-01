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
from scipy.special import erfcx,spence
from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.functions.log import log, exp
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

    def rankbound(self,Delta=1,function="sincsquared_fast"):
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
            sage: Z.rankbound(Delta=1,function="sincsquared_fast")
            2.0375000846
            sage: Z.rankbound(Delta=1,function="sincsquared")
            2.0375000846
            sage: Z.rankbound(Delta=1,function="gaussian")
            2.05689042503

        """

        # If Delta>6.95, then exp(2*pi*Delta)>sys.maxint, so we get overflow
        # when summing over the logarithmic derivative coefficients
        if Delta > 6.95:
            raise ValueError("Delta value too large; will result in overflow")

        if function=="sincsquared_fast":
            return self._rankbound_sincsquared_fast(Delta=Delta)
        elif function=="sincsquared":
            return self._rankbound_sincsquared(Delta=Delta)
        elif function=="gaussian":
            return self._rankbound_gaussian(Delta=Delta)
        else:
            raise ValueError("Input function not recognized.")

    def _rankbound_sincsquared(self,Delta=1):
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

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.rankbound`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z._rankbound_sincsquared(Delta=1)
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
            cn  = self.logarithmic_derivative_coefficient(n)
            if cn!=0:
                logn = log(RDF(n))
                y += cn*(t-logn)

        return 2*(u+w+y)/(t**2)

    def _rankbound_gaussian(self,Delta=1):
        r"""
        Bound from above the analytic rank of the form attached to self
        by computing
            '\sum_{\gamma} f(\Delta*\gamma),'
        where '\gamma' ranges over the imaginary parts of the zeros of 'L_E(s)'
        along the critical strip, and 'f(x) = exp(-x^2)'

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

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.rankbound`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z._rankbound_gaussian(Delta=1)
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
            cn  = self.logarithmic_derivative_coefficient(n)
            if cn != 0:
                logn = log(RDF(n))
                y += cn*exp(-(logn/(2*Delta))**2)

        # y is the truncation of an infinite sum, so we must add a value which
        # exceeds the max amount we could have left out.
        return RDF(u+w+y+0.1)/Deltasqrtpi

    def _rankbound_sincsquared_fast(self,Delta=1):
        """
        A faster, more intelligent implementation of self._rankbound_sincsquared().

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

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.rankbound_sincsquared`
            for the more general but slower version of this method.

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.rankbound`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z._rankbound_sincsquared_fast(Delta=1)
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

            sage: E = EllipticCurve([23,-100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.elliptic_curve()
            Elliptic Curve defined by y^2 = x^3 + 23*x - 100 over Rational Field

        """
        return self._E

    def logarithmic_derivative_coefficient(self,n):
        r"""
        Return the 'n'th derivative of the logarithmic derivative of the
        L-function attached to self, shifted so that the critical line
        lies on the imaginary axis. This is zero if 'n' is not a perfect
        prime power; when 'n=p^e' it is '-a_p^e*log(p)/p^e', where
        'a_p = p+1-\#{E(FF_p)}' is the trace of Frobenius at 'p'.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT:

        A real number at most '\frac{log(n)}{\sqrt{n}}' in magnitude

        EXAMPLES::

        sage: E = EllipticCurve('11a')
        sage: Z = LFunctionZeroSum(E)
        sage: for n in range(1,12): print(n,Z.logarithmic_derivative_coefficient(n))
        (1, 0)
        (2, 0.69314718056)
        (3, 0.366204096223)
        (4, 0.0)
        (5, -0.321887582487)
        (6, 0)
        (7, 0.555974328302)
        (8, -0.34657359028)
        (9, 0.610340160371)
        (10, 0)
        (11, -0.217990479345)

        """
        n = ZZ(n)
        if n==1:
            return ZZ(0)
        if not n.is_prime_power():
            return ZZ(0)

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
