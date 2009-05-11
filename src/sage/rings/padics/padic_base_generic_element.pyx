"""
p-Adic Base Generic Element.

A common superclass for features shared among all elements of Zp and
Qp (regardless of implementation).

AUTHORS::

    - David Roe
"""

from sage.rings.integer import Integer

cdef class pAdicBaseGenericElement(pAdicGenericElement):
    def __init__(self, parent):
        """
        Initialization

        EXAMPLES::

            sage: R = Zp(5); a = R(1287) #indirect doctest
        """
        self.prime_pow = <PowComputer_base>parent.prime_pow
        pAdicGenericElement.__init__(self, parent)

    cdef int _set_to_mpz(self, mpz_t dest) except -1:
        raise NotImplementedError

    cdef int _set_to_mpq(self, mpq_t dest) except -1:
        raise NotImplementedError

    def _pari_init_(self):
        """
        A string that when input into pari will yield the pari version
        of this element.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: pari(R(1777)) #indirect doctest
            2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
        """
        return self.lift().str() + " + O(" + self.parent().prime().str() + "^" + self.precision_absolute().str() + ")"

    def _integer_(self, Z=None):
        """
        Returns an integer congruent to this element modulo
        p^self.absolute_precision().

        EXAMPLES::

            sage: R = Zp(5); ZZ(R(-1)) #indirect doctest
            95367431640624
        """
        return Integer(self.lift())

    cdef int _cmp_units(left, pAdicGenericElement right):
        """
        Compares left and right, assuming both are units.

        TESTS::

            sage: R = Zp(5); a = R(5+5^7); b = R(5, 5); a == b #indirect doctest
            True
        """
        p = left.parent().prime()
        a = left.lift()
        b = right.lift()
        prec = min(left.precision_relative(), right.precision_relative())
        ppow = p**prec
        a %= ppow
        b %= ppow
        if a < b:
            return -1
        elif a == b:
            return 0
        else:
            return 1

    cpdef abs(self, prec=None):
        """
        Returns the p-adic absolute value of self.

        This is normalized so that the absolute value of p is 1/p.

        INPUT --
        prec - Integer.  The precision of the real field in which the answer is returned.  If None, returns a rational for absolutely unramified fields, or a real with 53 bits of precision if ramified.

        EXAMPLES:
        sage: a = Qp(5)(15); a.abs()
        1/5
        sage: a.abs(53)
        0.200000000000000
        """
        if prec is None:
            return self.prime_pow.prime**(-self.valuation())
        from sage.rings.real_mpfr import RealField
        return RealField(prec)(self.prime_pow.prime**(-self.valuation()))

    def exp(self):
        r"""
        Compute the p-adic exp of any element of $\Z_p$ where the
        series converges.

        EXAMPLES

        Borrowed from log.::

            sage: Z13 = Zp(13, 10, print_mode='series')
            sage: a = Z13(13 + 6*13**2 + 2*13**3 + 5*13**4 + 10*13**6 + 13**7 + 11*13**8 + 8*13**9).add_bigoh(10); a
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: a.exp()
            1 + 13 + O(13^10)
            sage: Q13 = Qp(13, 10, print_mode='series')
            sage: a = Q13(13 + 6*13**2 + 2*13**3 + 5*13**4 + 10*13**6 + 13**7 + 11*13**8 + 8*13**9).add_bigoh(10); a
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: a.exp()
            1 + 13 + O(13^10)

        The next few examples illustrate precision when computing
        $p$-adic exps.  First we create a field with \emph{default}
        precision 10.::

            sage: R = Zp(5,10, print_mode='series')
            sage: e = R(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*R.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        ::

            sage: K = Qp(5,10, print_mode='series')
            sage: e = K(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*K.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        TESTS

        Check that results are consistent over a range of precision::

            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_exp = (K(p)).exp()
            sage: for prec in range(2, max_prec):
            ...       ll = (K(p).add_bigoh(prec)).exp()
            ...       assert ll == full_exp
            ...       assert ll.precision_absolute() == prec
            sage: K = Qp(p, max_prec)
            sage: full_exp = (K(p)).exp()
            sage: for prec in range(2, max_prec):
            ...       ll = (K(p).add_bigoh(prec)).exp()
            ...       assert ll == full_exp
            ...       assert ll.precision_absolute() == prec

        AUTHORS::

            - Genya Zaytman (2007-02-15)

        """

        val = self.valuation()
        p = self.parent().prime()
        if val > 1 or (val > 0 and p != 2):

            prec = self.precision_absolute()

            # I believe this works
            max_term = ((p-1)*(prec-1))//((p-1)*val - 1) + 1

            # Need extra precision to take into account powers of p
            # in the denominators of the series. (Indeed, it's a
            # not-entirely-trivial fact that if x is given mod p^n, that
            # exp(x) is well-defined mod p^n !) .
            extra_prec = max_term//(p-1)

            from sage.rings.padics.factory import Zp
            working_ring = Zp(p, prec + extra_prec, type = 'capped-abs')
            x = working_ring(self.lift())
            term = ans = working_ring(Integer(1))
            for n in range(1, max_term):
                term *=x
                term = term // working_ring(Integer(n))
                ans += term
            # Note that it is the absolute precision that is respected by exp: even when p == 2?
            return self.parent()(ans).add_bigoh(prec)
        else:
            raise ValueError, "series doesn't converge"

    def log(self, branch = None):
        r"""
        Compute the p-adic logarithm of any unit in $\Z_p$.
        (See below for normalization.)

        The usual power series for log with values in the additive
        group of $\Z_p$ only converges for 1-units (units congruent to
        1 modulo p).  However, there is a unique extension of log to a
        homomorphism defined on all the units.  If u = a*v is a unit
        with v = 1 (mod p), then we define log(u) = log(v).  This is
        the correct extension because the units U of Z_p splits as a
        product U = V x <w>, where V is the subgroup of 1-units and w
        is a (p-1)st root of unity.  The <w> factor is torsion, so
        must go to 0 under any homomorphism to the torsion free group
        $(\Z_p, +)$.

        NOTES

        What some other systems do::

            * PARI:  Seems to define log the same way as we do.

            * MAGMA: Gives an error when unit is not a 1-unit.

        ALGORITHM

        Input: Some p-adic unit u.

        1. Check that the input p-adic number is really a unit (i.e.,
           valuation 0)

        2. Let $1-x = u^{p-1}$, which is a 1-unit.

        3. Use the series expansion

        ..math ::

            \log(1-x) = F(x) = -x - 1/2*x^2 - 1/3*x^3 - 1/4*x^4 - 1/5*x^5 - ...

        to compute the logarithm log(u**(p-1)).  Use enough terms so
        that terms added on are zero

        4. Then

        ..math ::

            \log(u) = log(u^{p-1})/(p-1) = F(1-u^{p-1})/(p-1).

        EXAMPLES::

            sage: Z13 = Zp(13, 10, print_mode='series')
            sage: a = Z13(14); a
            1 + 13 + O(13^10)

        Note that the relative precision decreases when we take log:
        it is the absolute precision that is preserved.::

            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: Q13 = Qp(13, 10, print_mode='series')
            sage: a = Q13(14); a
            1 + 13 + O(13^10)
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)

        The next few examples illustrate precision when computing
        $p$-adic logs.  First we create a field with \emph{default}
        precision 10.::

            sage: R = Zp(5,10, print_mode='series')
            sage: e = R(389); e
            4 + 2*5 + 3*5^3 + O(5^10)
            sage: e.log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: K = Qp(5,10, print_mode='series')
            sage: e = K(389); e
            4 + 2*5 + 3*5^3 + O(5^10)
            sage: e.log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)

        Check that results are consistent over a range of precision::

            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_log = (K(1 + p)).log()
            sage: for prec in range(2, max_prec):
            ...       ll = (K(1 + p).add_bigoh(prec)).log()
            ...       assert ll == full_log
            ...       assert ll.precision_absolute() == prec


        AUTHORS::

            - William Stein: initial version
            - David Harvey (2006-09-13): corrected subtle precision bug
              (need to take denominators into account! -- see trac \#53)
            - Genya Zaytman (2007-02-14): adapted to new p-adic class

        TODO:

            - Currently implemented as $O(N^2)$. This can be improved
              to soft-$O(N)$ using algorithm described by Dan
              Bernstein:
              http://cr.yp.to/lineartime/multapps-20041007.pdf
        """

        p = self.parent().prime()
        # Step 1 -- a unit?
        if self.is_unit() and self.unit_part().residue(1) == 1:
            # It's already a 1-unit, so just use the series
            # (base case of "induction")

            prec = self.precision_absolute()

            # Need extra precision to take into account powers of p
            # in the denominators of the series. (Indeed, it's a
            # not-entirely-trivial fact that if x is given mod p^n, that
            # log(x) is well-defined mod p^n !) Specifically:
            # we are only guaranteed that $x^j/j$ is zero mod $p^n$ if
            # j >= floor(log_p(j)) + n.
            extra_prec = 0
            while extra_prec < Integer(prec + extra_prec).exact_log(p):
                extra_prec += 1

            x = Integer(1) - self
            from sage.rings.padics.factory import Zp
            working_ring = Zp(p, prec + extra_prec, type = 'capped-abs', check=False)
            x = working_ring(x.lift())
            xpow = x
            ans = working_ring(Integer(0))
            for n in range(1, prec + extra_prec):
                ans -= xpow//working_ring(Integer(n))
                xpow *= x
            # Note that it is the absolute precision that is respected by log
            return self.parent()(ans.lift()).add_bigoh(prec)
        elif self.is_unit():
            return (self**Integer(p-1)).log() // Integer(p-1)
        elif not branch is None and self.parent().__contains__(branch):
            branch = self.parent()(branch)
            return self.unit_part().log() + branch*self.valuation()
        else:
            raise ValueError, "not a unit: specify a branch of the log map"

    def minimal_polynomial(self, name):
        """
        Returns a minimal polynomial of this p-adic element, i.e., x - self

        INPUT::

            - self -- a p-adic element
            - name -- string the name of the variable

        OUTPUT::

            polynomial -- a minimal polynomial of this p-adic element, i.e., x - self

        EXAMPLES::

            sage: Zp(5,5)(1/3).minimal_polynomial('x')
            (1 + O(5^5))*x + (3 + 5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5))
        """
        R = self.parent()[name]
        return R.gen() - R(self)

    def norm(self, ground=None):
        """
        Returns the norm of this p-adic element over the ground ring.

        NOTE!  This is not the p-adic absolute value.  This is a field
        theoretic norm down to a ground ring.  If you want the p-adic
        absolute value, use the abs() function instead.

        INPUT::

            - self -- a p-adic element
            - ground -- a subring of the ground ring (default: base ring)

        OUTPUT::

            element -- the norm of this p-adic element over the ground ring

        EXAMPLES::

            sage: Zp(5)(5).norm()
            5 + O(5^21)
        """
        if (ground is not None) and (ground != self.parent()):
            raise ValueError, "Ground Field not a subfield"
        else:
            return self

    def trace(self, ground=None):
        """
        Returns the trace of this p-adic element over the ground ring

        INPUT::

            - self -- a p-adic element
            - ground -- a subring of the ground ring (default: base ring)

        OUTPUT::

            element -- the trace of this p-adic element over the ground ring

        EXAMPLES::

            sage: Zp(5,5)(5).trace()
            5 + O(5^6)
        """
        if (ground is not None) and (ground != self.parent()):
            raise ValueError, "Ground ring not a subring"
        else:
            return self

    cdef int teichmuller_set_c(self, mpz_t value, mpz_t ppow) except -1:
        r"""
        Sets value to the integer between 0 and p^prec that is
        congruent to the Teichmuller lift of value to Z_p.  Does not
        affect self.

        INPUT::

            - value -- An mpz_t currently holding an approximation to
                       the Teichmuller representative (this
                       approximation can be any integer).  It will be
                       set to the actual Teichmuller lift

            - ppow -- An mpz_t holding the value p^prec, where prec is the
                      desired precision of the Teichmuller lift

        EXAMPLES::

            sage: R = Zp(next_prime(40000), 3); R.teichmuller(237) #indirect doctest
            237 + 7067*40009 + 5212*40009^2 + O(40009^3)
        """
        cdef mpz_t u, xnew
        if mpz_divisible_p(value, self.prime_pow.prime.value) != 0:
            mpz_set_ui(value, 0)
            return 0
        if mpz_sgn(value) < 0 or mpz_cmp(value, ppow) >= 0:
            mpz_mod(value, value, ppow)
        mpz_init(u)
        mpz_init(xnew)
        # u = 1 / Mod(1 - p, ppow)
        mpz_sub(u, ppow, self.prime_pow.prime.value)
        mpz_add_ui(u, u, 1)
        mpz_invert(u, u, ppow)
        # Consider x as Mod(self.value, ppow)
        # xnew = x + u*(x^p - x)
        mpz_powm(xnew, value, self.prime_pow.prime.value, ppow)
        mpz_sub(xnew, xnew, value)
        mpz_mul(xnew, xnew, u)
        mpz_add(xnew, xnew, value)
        mpz_mod(xnew, xnew, ppow)
        # while x != xnew:
        #     x = xnew
        #     xnew = x + u*(x^p - x)
        while mpz_cmp(value, xnew) != 0:
            mpz_set(value, xnew)
            mpz_powm(xnew, value, self.prime_pow.prime.value, ppow)
            mpz_sub(xnew, xnew, value)
            mpz_mul(xnew, xnew, u)
            mpz_add(xnew, xnew, value)
            mpz_mod(xnew, xnew, ppow)
        mpz_clear(u)
        mpz_clear(xnew)
