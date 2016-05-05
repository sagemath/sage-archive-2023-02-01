from types import MethodType

include "sage/libs/linkages/padics/fmpz_poly_unram.pxi"
include "sage/libs/linkages/padics/unram_shared.pxi"
include "CA_template.pxi"

cdef class PowComputer_(PowComputer_flint_unram):
    """
    A PowComputer for a capped-absolute unramified ring.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        """
        Initialization.

        EXAMPLES::

            sage: R.<a> = Zq(125)
            sage: type(R.prime_pow)
            <type 'sage.rings.padics.qadic_flint_CA.PowComputer_'>
            sage: R.prime_pow._prec_type
            'capped-abs'
        """
        self._prec_type = 'capped-abs'
        PowComputer_flint_unram.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)

cdef class qAdicCappedAbsoluteElement(CAElement):
    frobenius = MethodType(frobenius_unram, None, qAdicCappedAbsoluteElement)
    trace = MethodType(trace_unram, None, qAdicCappedAbsoluteElement)
    norm = MethodType(norm_unram, None, qAdicCappedAbsoluteElement)

    def matrix_mod_pn(self):
        """
        Returns the matrix of right multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for this
        extension field.  Thus the *rows* of this matrix give the
        images of each of the `x^i`.  The entries of the matrices are
        IntegerMod elements, defined modulo ``p^(self.absprec() / e)``.

        EXAMPLES::

            sage: R.<a> = ZqCA(5^5,5)
            sage: b = (5 + 15*a)^3
            sage: b.matrix_mod_pn()
            [ 125 1125  250  250    0]
            [   0  125 1125  250  250]
            [2375 2125  125 1125  250]
            [2375 1375 2125  125 1125]
            [2875 1000 1375 2125  125]

            sage: M = R(0,3).matrix_mod_pn(); M == 0
            True
            sage: M.base_ring()
            Ring of integers modulo 125
        """
        return cmatrix_mod_pn(self.value, self.absprec, 0, self.prime_pow)

    def _flint_rep(self, var='x'):
        """
        Replacement for _ntl_rep for use in printing and debugging.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, 4)
            sage: (1+a).inverse_of_unit()._flint_rep()
            41*x^2 + 40*x + 42
            sage: (1+a)*(41*a^2+40*a+42)
            1 + O(3^4)
        """
        return self.prime_pow._new_fmpz_poly(self.value, var)

    def _flint_rep_abs(self, var='x'):
        """
        Replacement for _ntl_rep_abs for use in printing and debugging.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, 4)
            sage: (3+3*a)._flint_rep_abs()
            (3*x + 3, 0)
        """
        return self._flint_rep(var), Integer(0)

    def __hash__(self):
        r"""
        Raise a ``TypeError`` since this element is not hashable
        (:trac:`11895`.)

        TESTS::

            sage: K.<a> = ZqCA(9)
            sage: hash(a)
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CA.qAdicCappedAbsoluteElement'

        """
        # Eventually, hashing will be disabled for all (non-fixed-mod) p-adic
        # elements (#11895), until then, we only to this for types which did
        # not support hashing before we switched some elements to FLINT
        raise TypeError("unhashable type: 'sage.rings.padics.qadic_flint_CA.qAdicCappedAbsoluteElement'")
