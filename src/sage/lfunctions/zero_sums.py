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
    def __init__(self):
        """
        Should never be called.
        """
        pass

    def rankbound(self,Delta=1,function="sincsquared"):
        """
        Bound from above the analytic rank of the form attached to self
        by computing
            sum_gamma f(Delta*gamma),
        where gamma ranges over the imaginary parts of the zeros of L_E(s)
        along the critical strip, and f(x) is an appropriate continuous
        L_2 function such that f(0)=1.

        As Delta increases this sum limits from above to the analytic rank
        of E, as f(0) = 1 is counted with multiplicity r, and the other terms
        go to 0.

        Current options for function are
         - "sincquared" -- f(x) = (sin(pi*x)/(pi*x))^2
         - "gaussian"   -- f(x) = exp(-x^2)
        """

        # If Delta>6.95, then exp(2*pi*Delta)>sys.maxint, so we get overflow
        # when summing over the logarithmic derivative coefficients
        if Delta > 6.95:
            raise ValueError("Delta value too large; will result in overflow")

        if function=="sincsquared":
            return self._rankbound_sincsquared(Delta=Delta)
        elif function=="gaussian":
            return self._rankbound_gaussian(Delta=Delta)
        else:
            raise ValueError("Input function not recognized.")

    def _rankbound_sincsquared(self,Delta=1):
        """
        Bound from above the analytic rank of the form attached to self
        by computing
            sum_gamma f(Delta*gamma),
        where gamma ranges over the imaginary parts of the zeros of L_E(s)
        along the critical strip, and f(x) = (sin(pi*x)/(pi*x))^2.
        As Delta increases this sum limits from above to the analytic rank
        of E, as f(0) = 1 is counted with multiplicity r, and the other terms
        go to 0.
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
            cn  = self.log_deriv_coeff(n)
            if cn!=0:
                logn = log(RDF(n))
                y += cn*(t-logn)

        return 2*(u+w+y)/(t**2)

    def _rankbound_gaussian(self,Delta=1):
        """
        Bound from above the analytic rank of the form attached to self
        by computing
            sum_gamma f(Delta*gamma),
        where gamma ranges over the imaginary parts of the zeros of L_E(s)
        along the critical strip, and f(x) = exp(-x^2).
        As Delta increases this sum limits from above to the analytic rank
        of E, as f(0) = 1 is counted with multiplicity r, and the other terms
        go to 0.
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
        while n < exp2piDelta:
            n += 1
            cn  = self.log_deriv_coeff(n)
            if cn != 0:
                logn = log(RDF(n))
                y += cn*exp(-(logn/(2*Delta))**2)

        return RDF(u+w+y+0.1)/Deltasqrtpi

    def _rankbound_sincsquared_fast(self,Delta=1):
        """
        A faster, more intelligent version of self._rankbound_sincsquared().

        Bound from above the analytic rank of the form attached to self
        by computing
            sum_gamma f(Delta*gamma),
        where gamma ranges over the imaginary parts of the zeros of L_E(s)
        along the critical strip, and f(x) = (sin(pi*x)/(pi*x))^2.
        As Delta increases this sum limits from above to the analytic rank
        of E, as f(0) = 1 is counted with multiplicity r, and the other terms
        go to 0.

        This will only produce correct output if self._E is given by its
        global minimal model, i.e. if self._E.is_minimal()==True.
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
        """

        self._E = E
        if N is not None:
            self._N = N
        else:
            self._N = self._E.conductor()
        # PARI minicurve for computing ap coefficients
        self._e = E.pari_mincurve()

        self._pi = RDF(pi)
        self._euler_gamma = RDF(euler_gamma)

    def __repr__(self):
        """
        Representation of self.
        """
        s = "Zero sum estimator for L-function attached to "
        return s+str(self._E)

    def E(self):
        """
        Return the elliptic curve associated with self.
        """
        return self._E

    def level(self):
        """
        Return the level of self, equal to the conductor of E.
        """
        return self._N

    def log_deriv_coeff(self,n):
        """
        Return the nth derivative of the logarithmic derivative of the
        L-function, shifted so that the critical line lies on the imaginary
        axis. This is zero if n is not a perfect prime power, and if n=p^e,
        then log_deriv_coeff(n) = -(a_p)^e*log(p)/p^e, where
        a_p = p+1-#{E(FF_p)} is the trace of Frobenius at p.
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
            #print(n,-ap*logn/n_float)
            return -ap*logn/n_float
        else:
            p,e = n.perfect_power()
            ap = self._E.ap(p)
            logp = log(RDF(p))
            if p.divides(self._N):
                return - ap**e*logp/n_float
            c = CDF(ap,(4*p-ap**2).sqrt())/2
            aq = (2*(c**e).real()).round()
            #print(n,-aq*logp/n_float)
            return -aq*logp/n_float

    def _sincsquared_summand(self,p,t,bound):
        """
        Return the summand for the summand corresponing to prime p
        of the sinc^2(Delta*gamma) sum over the zeros of the L_E(s),
        where t = 2*pi*Delta. The returned sum runs over all powers
        of p <= bound = floor(exp(t)).

        p must be a prime python int.
        """

        logp = log(RDF(p))
        ap = RDF(self._e.ellap(p))
        p_float = RDF(p)
        n = p_float**2
        if 2*logp > t:
            return (-ap*logp/p_float)*(t-logp)
        else:
            y = ap/p_float
            c = CDF(ap,(4*p_float-ap**2).sqrt())/2
            d = c
            logn = logp*2
            while logn < t:
                n = n*p_float
                d = d*c
                aq = RDF((2*d.real()).round())
                y += (aq/n)*(t-logn)

                logn += logp
            return -y*logp

def LFunctionZeroSum(X,*args,**kwds):
    """
    Constructor for the LFunctionZeroSum class.
    """

    # Here to avoid import recursion
    from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field

    if isinstance(X,EllipticCurve_rational_field):
        return LFunctionZeroSum_EllipticCurve(X,*args,**kwds)

    else:
        s = "Zero sum estimator class currently only "\
            +"implemented for elliptic curves over QQ."
        raise NotImplementedError(s)
