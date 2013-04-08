"""
`p`-Adic Base Generic Element

A common superclass for features shared among all elements of `\mathbb{Z}_p` and
`\mathbb{Q}_p` (regardless of implementation).

AUTHORS:

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

    cdef int _set_mpz_into(self, mpz_t dest) except -1:
        raise NotImplementedError

    cdef int _set_mpq_into(self, mpq_t dest) except -1:
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
        return "%s + O(%s^%s)" % (self.lift(), self.parent().prime(), self.precision_absolute())

    def _integer_(self, Z=None):
        """
        Returns an integer congruent to this element modulo
        ``p^self.absolute_precision()``.

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

    def minimal_polynomial(self, name):
        """
        Returns a minimal polynomial of this `p`-adic element, i.e., ``x - self``

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``name`` -- string: the name of the variable

        OUTPUT:

        - ``polynomial`` -- a minimal polynomial of this `p`-adic element,
          i.e., ``x - self``

        EXAMPLES::

            sage: Zp(5,5)(1/3).minimal_polynomial('x')
            (1 + O(5^5))*x + (3 + 5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5))
        """
        R = self.parent()[name]
        return R.gen() - R(self)

    def norm(self, ground=None):
        """
        Returns the norm of this `p`-adic element over the ground ring.

        NOTE!  This is not the `p`-adic absolute value.  This is a field
        theoretic norm down to a ground ring.  If you want the `p`-adic
        absolute value, use the ``abs()`` function instead.

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``ground`` -- a subring of the ground ring (default: base
          ring)

        OUTPUT:

        - element -- the norm of this `p`-adic element over the ground
          ring

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
        Returns the trace of this `p`-adic element over the ground ring

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``ground`` -- a subring of the ground ring (default: base
          ring)

        OUTPUT:

        - ``element`` -- the trace of this `p`-adic element over the
          ground ring

        EXAMPLES::

            sage: Zp(5,5)(5).trace()
            5 + O(5^6)
        """
        if (ground is not None) and (ground != self.parent()):
            raise ValueError, "Ground ring not a subring"
        else:
            return self

    def frobenius(self, arithmetic=True):
        """
        Applies a Frobenius automorphism to this element. This is the identity
        map, since this element lies in `\QQ_p`; it exists for compatibility
        with the
        :meth:`~sage.rings.padics.padic_ext_element.pAdicExtElement.frobenius`
        method of elements of extensions of `\QQ_p`.

        INPUT:

        - ``self`` -- a `p`-adic element
        - ``arithmetic`` -- whether to apply arithmetic Frobenius (as opposed
          to geometric Frobenius) -- ignored, since both are the identity map
          anyway.

        OUTPUT:

        - returns ``self``.

        EXAMPLES::

            sage: Qp(7)(2).frobenius()
            2 + O(7^20)
        """
        return self

    cdef int teichmuller_set_c(self, mpz_t value, mpz_t ppow) except -1:
        r"""
        Sets ``value`` to the integer between 0 and ``ppow`` that is
        congruent to the Teichmuller lift of value to `\mathbb{Z}_p`.
        Does not affect ``self``.

        INPUT:

        - ``value`` -- An ``mpz_t`` currently holding an approximation
          to the Teichmuller representative (this approximation can be
          any integer).  It will be set to the actual Teichmuller lift

        - ``ppow`` -- An ``mpz_t`` holding the value ``p^prec``, where
          ``prec`` is the desired precision of the Teichmuller lift

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
