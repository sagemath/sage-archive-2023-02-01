include "../../ext/stdsage.pxi"

from sage.misc.misc import ellipsis_iter
from sage.calculus.all import SR

from sage.rings.real_mpfi import RIF
from sage.rings.all import infinity

from sage.rings.arith import factorial

from itertools import imap, count, repeat

cdef class RigidAnalyticFunction_disc(SageObject):
    """
    A class for rigid analytic functions on an open rigid analytic disc (ie power series).  But these functions have additional convergence information attached to them, in addition to a mechanism to generate arbitrary coefficients.

    This class should eventually inherit from element and have parent the ring (or algebra) of rigid algebraic functions on the disc.  But I'm lazy right now.
    """
    def __init__(self, base_ring, log_of_radius, log_order, log_offset, coeffs):
        """
        Initializes a rigid analytic function defined on a disk.

        INPUT:
          base_ring -- p-adic field.  the base ring of this function.
          log_of_radius -- Integer or Rational.  The logarithm base p of the radius of convergence of self.
          log_order -- logarithmic order of self.  See below for definition.
          log_offset -- logarithmic offset of self.  See below for definition.
          coeffs -- an object giving the coefficients of self.  See below for specification.

        NOTES:
        If n = log_order, L = log_of_radius and C = log_offset, then the valuation of the ith coefficient
        (normalized so that the valuation of p is 1) must be bounded below by
        i*L - n*log_p(i) - C

        coeffs.__iter__() should return an iterator over the coefficients of self (an object that has a next method and an __iter__ method that returns self).
        coeffs[i] should return the coefficient of x^i,
        coeffs.denom_iter() returns an iterator over the denominators of the coefficients of self.
        coeffs.numer_iter() returns an iterator over the numerators of the coefficients of self.  The nth term of numer_iter divided by the nth term of denom_iter, cast into the base ring, should be the nth term of iter.
        coeffs.denom(i) returns the denominator of coeffs[i]
        coeffs.numer(i) returns the numerator of coeffs[i]
        coeffs.numer_all_one() returns True iff the numerators of all coefficients are 1.
        coeffs.denom_all_one() returns True iff the denominators of all coefficients are 1.
        """
        self._base_ring = base_ring
        self._log_of_radius = log_of_radius
        self._e = self._base_ring.e()
        self._p = self._base_ring.prime()
        self._RIF_p = RIF(self._p)
        self._SR_p = SR(self._p)
        self._meL = RIF(-self._log_of_radius * self._e)
        self._log_order = log_order
        self._en = RIF(self._log_order * self._e)
        self._log_offset = log_offset
        self._eC = RIF(self._log_offset * self._e)
        self._coeffs = coeffs

    def base_ring(self):
        """
        Returns the base ring of this function.

        Will eventually be replaced by element's method.
        """
        return self._base_ring

    def log_of_radius(self):
        """
        Returns log_p(r), where r is the radius of the disc on which this function converges.
        """
        return self._log_of_radius

    def radius(self):
        """
        Returns the radius of the disc on which this function converges.

        This radius may be larger than the radius of convergence of the parent.
        """
        return self._SR_p**(self._log_of_radius)

    def log_order(self):
        """
        Returns the logarithmic order of self.

        If n = self.log_order() and C = self.log_offset(), then the valuation of the ith coefficient
        (normalized so that the valuation of p is 1) is bounded below by
        i*L - n*log_p(i) - C
        """
        return self._log_order

    def log_offset(self):
        """
        Returns the logarithmic offset of self.

        If n = log_order, L = log_of_radius and C = log_offset, then the valuation of the ith coefficient
        (normalized so that the valuation of p is 1) must be bounded below by
        i*L - n*log_p(i) - C
        """
        return self._log_offset

    def coeffs(self):
        """
        Returns the object defining the coefficients of self.  This object should satisfy the following specifications.

        self.coeffs().__iter__() returns an iterator over the coefficients of self (an object that has a next method and iter returns self).  The coefficients should be elements of self.base_ring()
        self.coeffs()[i] returns the coefficient of x^i.  The coefficients should be elements of self.base_ring().
        self.coeffs().start_term() returns smallest i so that the coefficient of x^i is nonzero.
        self.coeffs().valuation(i) returns the valuation of self.coeffs()[i]
        self.coeffs().denom_iter() returns an iterator over the denominators of the coefficients of self.
        self.coeffs().numer_iter() returns an iterator over the numerators of the coefficients of self.  The nth term of numer_iter divided by the nth term of denom_iter, cast into the base ring, should be the nth term of iter.
        self.coeffs().denom(i) returns the denominator of self.coeffs()[i]
        self.coeffs().numer(i) returns the numerator of self.coeffs()[i]
        self.coeffs().numer_all_one() returns True iff the numerators of all coefficients are 1 after the start_term.
        self.coeffs().denom_all_one() returns True iff the denominators of all coefficients are 1 after the start term.
        """
        return self._coeffs

    cpdef Integer coeff_val_bound(self, i):
        """
        Returns a non-strict lower bound on the valuation of the ith coefficient of self.
        """
        return self.coeff_val_bound_RIF(i).lower().floor()

    cpdef RealIntervalFieldElement coeff_val_bound_RIF(self, _i):
        """
        Returns a non-strict lower bound on the valuation of the ith coefficient of self, as an RIF.
        """
        cdef RealIntervalFieldElement i
        if PY_TYPE_CHECK(_i, RealIntervalFieldElement):
            i = <RealIntervalFieldElement>_i
        else:
            i = RIF(_i)
        return -i * self._meL - self._en * i.log(self._RIF_p) - self._eC

    cpdef Integer term_val_bound(self, i, val):
        """
        Returns a lower bound on the valuation of a_i * x^i when x has valuation i.
        """
        return self.term_val_bound_RIF(i, val).lower().floor()

    cpdef RealIntervalFieldElement term_val_bound_RIF(self, _i, val):
        """
        Returns a lower bound on the valuation of a_i * x^i when x has valuation i, as an RIF.

        If val is infinity, returns RIF(0)
        """
        if val is infinity:
            return RIF(0)
        return self.coeff_val_bound_RIF(_i) + RIF(val * _i)

    def needed_indices(self, absprec, val = None, empty_list = False):
        """
        Returns the indices needed to compute the value of this function to absolute precision absprec if the input has valuation val.

        INPUT:
          self -- a RigidAnalyticFunction_disc
          absprec -- the absolute precision desired of the answer.
          val -- the valuation of the element which is being plugged in (normalized so that the valuation of the uniformizer of the base ring is 1) [default None].  Defaults to the minimum valuation that's strictly greater than -self.log_of_radius() * self.base_ring().e()
          empty_list -- boolean (default False).  Whether to use self.coeffs() to decrease max_contiguous and put outliers into L.  If True, does not do so.
        OUTPUT:
          max_contiguous -- non-negative integer.  Use the iterator for indices 0 <= i < max_contiguous
          L -- list of non-negative integers.  The additional terms needed that affect the return value below absolute precision absprec.
        """
        if val is None:
            val = (self._meL + 1).lower().floor()
        elif val <= self._meL:
            raise ValueError, "val must be greater than %s"%self._meL
        elif val is infinity:
            return 1, []
        # check the following
        cdef Integer ans
        # There's surely a faster way to find this.
        ans = ((absprec - self._log_offset) / (val - self._meL)).upper().ceil()
        while (ans * (val - self._meL) - self._en * RIF(ans).log(self._RIF_p) - self._eC).lower().floor() < absprec:
            ans += 1
        orig = ans
        L = []
        if not empty_list:
            ans_reset = False
            i = Integer(self._coeffs.start_term())
            for a in self._coeffs:
                #print "a.valuation = %s"%a.valuation()
                #print "i = %s, val = %s"%(i, a.valuation() + i *val)
                if a.valuation() + i * val >= absprec:
                    if not ans_reset:
                        ans_reset = True
                        ans = i
                else:
                    if ans_reset:
                        L.append(i)
                if i >= orig:
                    break
                i += 1
        return ans, L

    cpdef output_valuation(self, input_valuation):
        """
        Returns the valuation of this function when evaluated on something of valuation input_valuation.

        Neglects possible cancellation (so this is a lower bound on the valuation of the answer) and with valuation normalized so that the valuation of the uniformizer is 1.
        """
        ## We have that if this function is given by the power series $\sum_{i = 0}^{\infty} a_i x^i$, then the valuation of a_i x^i
        return self.valuation_low_point_pair(input_valuation)[1]

    cpdef output_valuation_fast(self, input_valuation):
        """
        Returns a lower bound on the valuation of this function when evaluated on something of valuation input_valuation.
        """
        if input_valuation is infinity:
            if self._coeffs.start_term() == 0:
                return self._coeffs.valuation(0)
            else:
                return infinity
        cdef RealIntervalFieldElement minimum = self.valuation_low_point_fast(input_valuation)
        return min([self._coeffs.valuation(i) + i * input_valuation for i in range(minimum.lower().floor(), minimum.upper().ceil() + 1)])

    cpdef valuation_low_point(self, input_valuation):
        """
        Returns an index i that minimizes the valuation of the term a_i * x^i given that x had valuation input_valuation.
        """
        return self.valuation_low_point_pair(input_valuation)[0]

    cpdef valuation_low_point_pair(self, input_valuation):
        """
        Returns a pair giving an index i and valuation v so that v is minimal among valuations of terms a_i * x^i where x has valuation input_valuation.
        """
        if input_valuation is infinity:
            if self._coeffs.start_term() == 0:
                return 0, self._coeffs.valuation(0)
            else:
                return 0, infinity
        cdef Integer initial = self.valuation_low_point_fast(input_valuation).lower().floor()
        start = self._coeffs.start_term()
        if initial < start:
            initial = start
        cdef Integer candidate = initial
        min_val_found = self._coeffs.valuation(candidate) + candidate * input_valuation
        cdef Integer i = initial - 1
        while i >= start and self.term_val_bound_RIF(i, input_valuation) < min_val_found:
            test_val_found = self._coeffs.valuation(i) + i * input_valuation
            if test_val_found < min_val_found:
                min_val_found = test_val_found
                candidate = i
            i -= 1
        i = initial + 1
        while self.term_val_bound_RIF(i, input_valuation) < min_val_found:
            test_val_found = self._coeffs.valuation(i) + i * input_valuation
            if test_val_found < min_val_found:
                min_val_found = test_val_found
                candidate = i
            i += 1
        return candidate, min_val_found

    cpdef RealIntervalFieldElement valuation_low_point_fast(self, input_valuation):
        """
        Returns an index i that is within 1 of the one that minimizes the function i * (input_valuation + e * log_p(r)) - e * n * log_p(i) - e * C

        If input_valuation is infinity, returns RIF(0.5)
        """
        # Differentiating yields input_valuation + e * log_p(r) = e * n / (ln(p) * i)
        # so i = e * n / (ln(p) * (input_valuation + e * log_p(r)))
        if input_valuation is infinity:
            return RIF(0.5)
        return self._en / (self._RIF_p.log() * (RIF(input_valuation) - self._meL))

    def output_absprec(self, input_absprec, val = None):
        """
        """
        if val is None:
            val = (self._meL + 1).lower().floor()
        elif val <= self._meL:
            raise ValueError, "val must be greater than %s"%self._meL
        elif val > input_absprec:
            raise ValueError, "val must be less than absprec"
        elif val == input_absprec:
            return self.output_valuation(val)
        input_relprec = input_absprec - val
        cdef Integer initial = self.valuation_low_point_fast(val).lower().floor()
        start = self._coeffs.start_term()
        if self._p.divides(initial):
            initial += 1
        if initial < start:
            initial = start
        min_aprec_found = self._coeffs.valuation(candidate) + candidate * val + input_relprec
        cdef Integer candidate = initial
        cdef Integer i = initial - 1
        while i >= start and self.term_val_bound_RIF(i, val) + input_relprec < min_aprec_found:
            ival = i.valuation(self._p)
            test_relprec = input_relprec
            while ival > 0 and test_relprec < self._e / (self._p - 1):
                test_relprec *= self._p
                ival -= 1
            if ival > 0:
                test_relprec += self._e * ival
            test_aprec_found = self._coeffs.valuation(i) + i * val + test_relprec
            if test_aprec_found < min_aprec_found:
                min_aprec_found = test_aprec_found
                candidate = i
            i -= 1
        i = initial + 1
        while self.term_val_bound_RIF(i, val) + input_relprec < min_aprec_found:
            ival = i.valuation(self._p)
            test_relprec = input_relprec
            while ival > 0 and test_relprec < self._e / (self._p - 1):
                test_relprec *= self._p
                ival -= 1
            if ival > 0:
                test_relprec += self._e * ival
            test_aprec_found = self._coeffs.valuation(i) + i * val + test_relprec
            if test_aprec_found < min_aprec_found:
                min_aprec_found = test_aprec_found
                candidate = i
            i += 1
        return candidate, min_aprec_found

cdef class RigidAnalyticFunctionCoeffs(SageObject):
    def __init__(self, base_ring):
        self._base_ring = base_ring

    def __iter__(self):
        return imap(self.__getitem__, count(self.start_term()))

    def __getitem__(self, i):
        raise NotImplementedError

    def start_term(self):
        return Integer(0)

    def valuation(self, i):
        return self[i].valuation()

    def denom_iter(self):
        return repeat(self._base_ring(1))

    def numer_iter(self):
        return self.__iter__()

    def denom(self, i):
        return self._base_ring(1)

    def numer(self, i):
        return self[i]

    def numer_all_one(self):
        return False

    def denom_all_one(self):
        return True

cdef class ExpCoeffs(RigidAnalyticFunctionCoeffs):
    def __getitem__(self, i):
        return ~self._base_ring(factorial(i))

cdef class LogCoeffs(RigidAnalyticFunctionCoeffs):
    def __getitem__(self, i):
        return ~self._base_ring(i)

    def start_term(self):
        return Integer(1)

cdef class RigidAnalyticFunction_exp(RigidAnalyticFunction_disc):
    def __init__(self, base_ring):
        RigidAnalyticFunction_disc.__init__(self, base_ring, -~(base_ring.prime() - 1), 0, 0, ExpCoeffs(base_ring))


cdef class RigidAnalyticFunction_log(RigidAnalyticFunction_disc):
    def __init__(self, base_ring):
        RigidAnalyticFunction_disc.__init__(self, base_ring, 0, 1, RIF(0.001), LogCoeffs(base_ring))

    cpdef Integer coeff_val_bound(self, _i):
        cdef Integer i
        if PY_TYPE_CHECK(_i, Integer):
            i = <Integer> _i
        else:
            i = Integer(i)
        return i.valuation(self._p) * self._e

    cpdef Integer term_val_bound(self, _i, val):
        return self.coeff_val_bound(_i) + val * _i

    cpdef output_valuation(self, input_valuation):
        pass
