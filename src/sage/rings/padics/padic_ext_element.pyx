"""
p-Adic Extension Element.

A common superclass for all elements of extension rings and field of Zp and Qp.

AUTHORS:

- David Roe
"""

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.integer import Integer
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p

cdef class pAdicExtElement(pAdicGenericElement):
    cdef int _set_from_list(self, L) except -1:
        """
        Sets self from a list.

        The list should either be uniform in type, or all of the entries should be coercible to integers.
        If any of the entries in L is a list, L will be cast to a ZZ_pEX

        INPUT:
        L -- a list.
        """
        raise NotImplementedError

    cdef int _set_from_list_rel(self, L, long relprec) except -1:
        raise NotImplementedError

    cdef int _set_from_list_abs(self, L, long absprec) except -1:
        raise NotImplementedError

    cdef int _set_from_list_both(self, L, long absprec, long relprec) except -1:
        raise NotImplementedError

    cdef int _set_from_ZZX(self, ZZX_c poly) except -1:
        """
        Sets from a ZZX_c, choosing how to handle based on the
        precision type of self.parent().

        Fixed modulus elements should override this function.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZX_rel(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZX_abs(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZX should have been overridden"

    cdef int _set_from_ZZX_rel(self, ZZX_c poly, long relprec) except -1:
        """
        Set from a ZZX_c with bounded relative precision.

        Capped relative elements should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZX_both(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZX_abs(self, ZZX_c poly, long absprec) except -1:
        """
        Set from a ZZX_c with bounded absolute precision.

        Capped absolute elements should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZX_both(poly, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, long relprec) except -1:
        """
        Set from a ZZX_c with both absolute and relative precisions bounded.

        This function should be overridden for both capped absolute
        and capped relative elements.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZX(poly)
        else:
            raise RuntimeError, "_set_from_ZZX_both should have been overridden"

    cdef int _set_from_ZZ_pX(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx) except -1:
        """
        Sets self from a ZZ_pX defined with context ctx.

        This function should be overridden for fixed modulus elements.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pX_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pX_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZ_pX should have been overridden"

    cdef int _set_from_ZZ_pX_rel(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pX_c with bounded relative precision.

        Capped relative rings should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZ_pX_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pX_abs(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pX_c with bounded absolute precision.

        Capped absolute rings should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZ_pX_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pX_both(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec, long relprec) except -1:
        """
        Set from a ZZ_pX_c with both absolute and relative precision bounded.

        This function should be overridden by both capped absolute and capped relative elements.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pX(poly, ctx)
        else:
            raise RuntimeError, "_set_from_ZZ_pX_both should have been overridden"

    cdef int _set_from_ZZ_pE(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx) except -1:
        """
        Set from a ZZ_pE_c.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pE_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pE_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZ_pE should have been overridden"

    cdef int _set_from_ZZ_pE_rel(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pE_c with bounded relative precision.

        Capped relative rings should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZ_pE_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pE_abs(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pE_c with bounded absolute precision.

        Capped absolute elements should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZ_pE_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pE_both(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1:
        """
        Sets from a ZZ_pE_c with both absolute and relative precision bounded.

        Capped absolute and capped relative elements should override
        this function.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pE(poly, ctx)
        else:
            raise RuntimeError, "_set_from_ZZ_pE_both should have been overridden"

    cdef int _set_from_ZZ_pEX(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx) except -1:
        """
        Sets self from a ZZ_pEX_c.

        Fixed modulus elements should override this function.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pEX_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pEX_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZ_pEX should have been overridden"

    cdef int _set_from_ZZ_pEX_rel(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pEX_c with bounded relative precision.

        Capped relative elements should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZ_pEX_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pEX_abs(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pEX_c with bounded absolute precision.

        Capped absolute elements should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZ_pEX_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pEX_both(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1:
        """
        Sets from a ZZ_pEX_c with both absolute and relative precision bounded.

        Capped absolute and capped relative elements should override
        this function.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pEX(poly, ctx)
        else:
            raise RuntimeError, "_set_from_ZZ_pEX_both should have been overridden"

    cdef long _check_ZZ_pContext(self, ntl_ZZ_pContext_class ctx) except -1:
        raise NotImplementedError

    cdef long _check_ZZ_pEContext(self, ntl_ZZ_pEContext_class ctx) except -1:
        raise NotImplementedError

    cdef ext_p_list(self, bint pos):
        raise NotImplementedError

    cdef ext_p_list_precs(self, bint pos, long prec):
        raise NotImplementedError

    def _const_term_test(self):
        """
        Returns the constant term of a polynomial representing self.

        This function is mainly for troubleshooting, and the meaning
        of the return value will depend on whether self is capped
        relative or otherwise.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566)
            sage: a._const_term_test()
            566
        """
        cdef ntl_ZZ_p ans = ntl_ZZ_p(modulus=self.parent().prime())
        ans.x = self._const_term()
        return ans

    cdef ZZ_p_c _const_term(self):
        raise NotImplementedError

    def _ext_p_list(self, pos):
        """
        Returns a list of integers (in the eisenstein case) or a list
        of lists of integers (in the unramified case).  self can be
        reconstructed as a sum of elements of the list times powers of
        the uniformiser (in the eisenstein case), or as a sum of
        powers of the p times polynomials in the generator (in the
        unramified case).

        Note that zeros are truncated from the returned list, so you
        must use the valuation() function to completely recover self.

        INPUTS::

            - pos -- bint.  If True, all integers will be in the range [0,p-1],
              otherwise they will be in the range [(1-p)/2, p/2].

        OUTPUT::

            - L -- A list of integers or list of lists giving the
              series expansion of self.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W(775, 19); y
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: y._ext_p_list(True)
            [1, 0, 4, 0, 2, 1, 2, 4, 1]
            sage: y._ext_p_list(False)
            [1, 0, -1, 0, 2, 1, 2, 0, 1]
        """
        return self.ext_p_list(pos)

    cpdef abs(self, prec=None):
        """
        Returns the p-adic absolute value of self.

        This is normalized so that the absolute value of p is 1/p.

        INPUT --
        prec - Integer.  The precision of the real field in which the answer is returned.  If None, returns a rational for absolutely unramified fields, or a real with 53 bits of precision if ramified.

        EXAMPLES:
        sage: R = Zp(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: w.abs()
        0.724779663677696
        """
        if prec is None:
            if self.parent().e() == 1:
                return self.parent().prime()**(-self.valuation())
            else:
                prec = 53
        from sage.rings.real_mpfr import RealField
        return RealField(prec)(self.parent().prime())**(-self.ordp())
