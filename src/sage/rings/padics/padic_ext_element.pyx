from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.integer import Integer


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
        if self.parent().is_capped_relative():
            self._set_from_ZZX_rel(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZX_abs(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZX should have been overridden"

    cdef int _set_from_ZZX_rel(self, ZZX_c poly, long relprec) except -1:
        """
        Set from a ZZX_c with bounded relative precision.

        Capped relative rings should override this function, so the default implementation is
        for capped absolute.
        """
        self._set_from_ZZX_both(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZX_abs(self, ZZX_c poly, long absprec) except -1:
        """
        Set from a ZZX_c with bounded absolute precision.

        Capped absolute rings should override this function, so the default implementation is
        for capped relative.
        """
        self._set_from_ZZX_both(poly, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, long relprec) except -1:
        if self.parent().is_fixed_mod():
            self._set_from_ZZX(poly)
        else:
            raise RuntimeError, "_set_from_ZZX_both should have been overridden"

    cdef int _set_from_ZZ_pX(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx) except -1:
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pX_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pX_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZ_pX should have been overridden"

    cdef int _set_from_ZZ_pX_rel(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pX_c with bounded relative precision.

        Capped relative rings should override this function, so the default implementation is
        for capped absolute.
        """
        self._set_from_ZZ_pX_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pX_abs(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pX_c with bounded absolute precision.

        Capped absolute rings should override this function, so the default implementation is
        for capped relative.
        """
        self._set_from_ZZ_pX_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pX_both(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec, long relprec) except -1:
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pX(poly, ctx)
        else:
            raise RuntimeError, "_set_from_ZZ_pX_both should have been overridden"

    cdef int _set_from_ZZ_pE(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx) except -1:
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pE_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pE_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZ_pE should have been overridden"

    cdef int _set_from_ZZ_pE_rel(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pE_c with bounded relative precision.

        Capped relative rings should override this function, so the default implementation is
        for capped absolute.
        """
        self._set_from_ZZ_pE_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pE_abs(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pE_c with bounded absolute precision.

        Capped absolute rings should override this function, so the default implementation is
        for capped relative.
        """
        self._set_from_ZZ_pE_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pE_both(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1:
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pE(poly, ctx)
        else:
            raise RuntimeError, "_set_from_ZZ_pE_both should have been overridden"

    cdef int _set_from_ZZ_pEX(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx) except -1:
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pEX_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pEX_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError, "_set_from_ZZ_pEX should have been overridden"

    cdef int _set_from_ZZ_pEX_rel(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pEX_c with bounded relative precision.

        Capped relative rings should override this function, so the default implementation is
        for capped absolute.
        """
        self._set_from_ZZ_pEX_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pEX_abs(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pEX_c with bounded absolute precision.

        Capped absolute rings should override this function, so the default implementation is
        for capped relative.
        """
        self._set_from_ZZ_pEX_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pEX_both(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1:
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

    cdef ZZ_p_c _const_term(self):
        raise NotImplementedError

    def _ext_p_list(self, pos):
        return self.ext_p_list(pos)
