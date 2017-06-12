#from types import MethodType

include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "FM_template.pxi"

cdef class RelativeRamifiedFixedModElement(FMElement):
    #frobenius = MethodType(frobenius_unram, None, qAdicFixedModElement)
    #trace = MethodType(trace_unram, None, qAdicFixedModElement)
    #norm = MethodType(norm_unram, None, qAdicFixedModElement)

    def __cinit__(self, parent=None, x=None, absprec=infinity, relprec=infinity):
        # It's not possible to set self.value in cconstruct (because of the calling syntax)
        # so we do it here.
        cdef type t
        if parent is not None: # This will break the pickling function
            t = type((<PowComputer_?>parent.prime_pow).modulus)
            self.value = t.__new__(t)
            #self.value = celement.__new__(celement)

    cdef FMElement _new_c(self):
        """
        Creates a new element with the same basic info.

        TESTS::

            sage: R = ZpFM(5); R(6) * R(7) #indirect doctest
            2 + 3*5 + 5^2 + O(5^20)
        """
        cdef type t = type(self)
        cdef type polyt = type(self.prime_pow.modulus)
        cdef FMElement ans = t.__new__(t)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        ans.value = polyt.__new__(polyt)
        #ans.value = celement.__new__(celement)
        cconstruct(ans.value, ans.prime_pow)
        return ans

    def _poly_rep(self):
        return self.value

    #def matrix_mod_pn(self):
    #    """
    #    Returns the matrix of right multiplication by the element on
    #    the power basis `1, x, x^2, \ldots, x^{d-1}` for this
    #    extension field.  Thus the *rows* of this matrix give the
    #    images of each of the `x^i`.  The entries of the matrices are
    #    IntegerMod elements, defined modulo ``p^(self.absprec() / e)``.

    #    EXAMPLES::

    #        sage: R.<a> = ZqFM(5^5,5)
    #        sage: b = (5 + 15*a)^3
    #        sage: b.matrix_mod_pn()
    #        [ 125 1125  250  250    0]
    #        [   0  125 1125  250  250]
    #        [2375 2125  125 1125  250]
    #        [2375 1375 2125  125 1125]
    #        [2875 1000 1375 2125  125]

    #        sage: M = R(0,3).matrix_mod_pn(); M == 0
    #        True
    #        sage: M.base_ring()
    #        Ring of integers modulo 3125
    #    """
    #    return cmatrix_mod_pn(self.value, self.prime_pow.prec_cap, 0, self.prime_pow)

    #def _flint_rep(self, var='x'):
    #    """
    #    Replacement for _ntl_rep for use in printing and debugging.

    #    EXAMPLES::

    #        sage: R.<a> = ZqFM(27, 4)
    #        sage: (1+a).inverse_of_unit()._flint_rep()
    #        41*x^2 + 40*x + 42
    #        sage: (1+a)*(41*a^2+40*a+42)
    #        1 + O(3^4)
    #    """
    #    return self.prime_pow._new_fmpz_poly(self.value, var)

    #def _flint_rep_abs(self, var='x'):
    #    """
    #    Replacement for _ntl_rep_abs for use in printing and debugging.

    #    EXAMPLES::

    #        sage: R.<a> = ZqFM(27, 4)
    #        sage: (3+3*a)._flint_rep_abs()
    #        (3*x + 3, 0)
    #    """
    #    return self._flint_rep(var), Integer(0)
