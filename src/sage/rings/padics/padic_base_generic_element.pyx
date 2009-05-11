

cdef class pAdicBaseGenericElement(pAdicGenericElement):
    def __init__(self, parent):
        self.prime_pow = <PowComputer_base>parent.prime_pow
        pAdicGenericElement.__init__(self, parent)

    cdef int _set_to_mpz(self, mpz_t dest) except -1:
        raise NotImplementedError

    cdef int _set_to_mpq(self, mpq_t dest) except -1:
        raise NotImplementedError

    cdef int _cmp_units(left, pAdicGenericElement right):
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

    cdef int teichmuller_set_c(self, mpz_t value, mpz_t ppow) except -1:
        r"""
        Sets value to the integer between 0 and p^prec that is congruent to the Teichmuller lift of value to Z_p.
        Does not affect self.

        INPUT:
        value -- An mpz_t currently holding an approximation to the
                 Teichmuller representative (this approximation can
                 be any integer).  It will be set to the actual
                 Teichmuller lift
        ppow -- An mpz_t holding the value p^prec, where prec is the desired precision of the Teichmuller lift
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
