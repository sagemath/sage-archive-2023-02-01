"""
`p`-Adic Base Generic Element

A common superclass for features shared among all elements of `\mathbb{Z}_p` and
`\mathbb{Q}_p` (regardless of implementation).

AUTHORS:

- David Roe
"""

from sage.rings.integer import Integer
from sage.misc.all import prod


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
        return "%s + O(%s^%s)" % (self.lift(), self.parent().prime(),
                                  self.precision_absolute())

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

        .. WARNING::

            This is not the `p`-adic absolute value. This is a field
            theoretic norm down to a ground ring. If you want the `p`-adic
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
            raise ValueError("Ground Field not a subfield")
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
            raise ValueError("Ground ring not a subring")
        else:
            return self

    def frobenius(self, arithmetic=True):
        """
        Applies a Frobenius automorphism to this element.

        This is the identity map, since this element lies in `\QQ_p`; it
        exists for compatibility with the
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

    def dwork_expansion(self, bd=20):
        r"""
        Return the value of a function defined by Dwork, used to compute
        the `p`-adic Gamma function.

        INPUT:

        - ``bd`` is a bound for precision, defaults to 20

        OUTPUT:

        - a ``p``-adic integer

        .. NOTE::

            This is based on GP code written by Fernando Rodriguez
            Villegas (http://www.ma.utexas.edu/cnt/cnt-frames.html).
            William Stein sped it up for GP
            (http://sage.math.washington.edu/home/wstein/www/home/wbhart/pari-2.4.2.alpha/src/basemath/trans2.c).
            The output is a `p`-adic integer from Dwork's expansion,
            used to compute the `p`-adic gamma function as in [RV]_
            section 6.2.

        REFERENCES:

        .. [RV] Rodriguez Villegas, Fernando. Experimental Number Theory.
           Oxford Graduate Texts in Mathematics 13, 2007.

        EXAMPLES::

            sage: R = Zp(17)
            sage: x = R(5+3*17+13*17^2+6*17^3+12*17^5+10*17^(14)+5*17^(17)+O(17^(19)))
            sage: x.dwork_expansion(18)
            16 + 7*17 + 11*17^2 + 4*17^3 + 8*17^4 + 10*17^5 + 11*17^6 + 6*17^7 + 17^8 + 8*17^10 + 13*17^11 + 9*17^12 + 15*17^13  + 2*17^14 + 6*17^15 + 7*17^16 + 6*17^17 + O(17^18)

            sage: R = Zp(5)
            sage: x = R(3*5^2+4*5^3+1*5^4+2*5^5+1*5^(10)+O(5^(20)))
            sage: x.dwork_expansion()
            4 + 4*5 + 4*5^2 + 4*5^3 + 2*5^4 + 4*5^5 + 5^7 + 3*5^9 + 4*5^10 + 3*5^11 + 5^13 + 4*5^14 + 2*5^15 + 2*5^16 + 2*5^17 + 3*5^18 + O(5^20)
        """
        R = self.parent()
        p = R.prime()
        s = R.one().add_bigoh(bd)
        t = s
        u = [s]
        for j in range(1, p):
            u.append(u[j-1] / j)
        for k in range(1, bd):
            u[0] = ((u[-1] + u[0]) / k) >> 1
            for j in range(1, p):
                u[j] = (u[j-1] + u[j]) / (j + k * p )
            t *= (self + k - 1)
            s += t * (u[0] << k)
        return R(-s)

    def padic_gamma(self, flag='pari'):
        r"""
        Return the value of the `p`-adic Gamma function.

        INPUT:

        - flag can be set to 'pari' to call the pari function, or 'sage' to call the function implemented in sage.
          set to 'pari' by default, since pari is about 10 times faster than sage.

        OUTPUT:

        - a `p`-adic integer

        .. NOTE::

            This is based on GP code written by Fernando Rodriguez
            Villegas (http://www.ma.utexas.edu/cnt/cnt-frames.html).
            William Stein sped it up for GP
            (http://sage.math.washington.edu/home/wstein/www/home/wbhart/pari-2.4.2.alpha/src/basemath/trans2.c).
            The 'sage' version uses dwork_expansion() to compute the
            `p`-adic gamma function of self as in [RV]_ section 6.2.

        EXAMPLES:

        This example illustrates ``x.padic_gamma()`` for `x` a `p`-adic unit.

        ::

            sage: R = Zp(7)
            sage: x = R(2+3*7^2+4*7^3+O(7^20))
            sage: x.padic_gamma('pari')
            1 + 2*7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^8 + 7^9 + 4*7^10 + 3*7^12 + 7^13 + 5*7^14 + 3*7^15 + 2*7^16 + 2*7^17 + 5*7^18 + 4*7^19 + O(7^20)
            sage: x.padic_gamma('sage')
            1 + 2*7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^8 + 7^9 + 4*7^10 + 3*7^12 + 7^13 + 5*7^14 + 3*7^15 + 2*7^16 + 2*7^17 + 5*7^18 + 4*7^19 + O(7^20)
            sage: x.padic_gamma('pari') == x.padic_gamma('sage')
            True

        Now ``x.padic_gamma()`` for `x` a `p`-adic integer but not a unit.

        ::

            sage: R = Zp(17)
            sage: x = R(17+17^2+3*17^3+12*17^8+O(17^13))
            sage: x.padic_gamma('pari')
            1 + 12*17 + 13*17^2 + 13*17^3 + 10*17^4 + 7*17^5 + 16*17^7 + 13*17^9 + 4*17^10 + 9*17^11 + 17^12 + O(17^13)
            sage: x.padic_gamma('sage')
            1 + 12*17 + 13*17^2 + 13*17^3 + 10*17^4 + 7*17^5 + 16*17^7 + 13*17^9 + 4*17^10 + 9*17^11 + 17^12 + O(17^13)
            sage: x.padic_gamma('pari') == x.padic_gamma('sage')
            True

        Finally, this function is not defined if `x` is not a `p`-adic integer.

        ::

            sage: K = Qp(7)
            sage: x = K(7^-5 + 2*7^-4 + 5*7^-3 + 2*7^-2 + 3*7^-1 + 3 + 3*7 + 7^3 + 4*7^4 + 5*7^5 + 6*7^8 + 3*7^9 + 6*7^10 + 5*7^11 + 6*7^12 + 3*7^13 + 5*7^14 + O(7^15))
            sage: x.padic_gamma()
            Traceback (most recent call last):
            ...
            ValueError: The p-adic gamma function only works on elements of Zp
        """
        if self.valuation() < 0:
            raise ValueError('The p-adic gamma function only works '
                             'on elements of Zp')
        else:
            if flag == 'pari':
                return self._pari_().gamma()

            elif flag == 'sage':
                p = self.parent().prime()
                n = self.precision_absolute()
                bd = n + 2*n//p
                if self.is_padic_unit():
                    k = Integer(self.residue())  # leading coefficient of self, non-zero
                    x = (self-k) >> 1
                    return (-1)**(k+1)*x.dwork_expansion(bd)*prod(j + (x << 1) for j in range(1, k))
                else:
                    return -(self >> 1).dwork_expansion(bd)
