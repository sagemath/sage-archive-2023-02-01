include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "FP_template.pxi"

cdef class RelativeRamifiedFloatingPointElement(FPElement):
    def _poly_rep(self):
        """
        Return the underlying polynomial representation of the element
        (which is used for computations).

        For debugging and printing purpose.

        EXAMPLES::

            sage: K.<a> = QqFP(125)
            sage: S.<x> = PolynomialRing(K)
            sage: W.<w> = K.extension(x^3 - 25*x^2 - 5*a*x + 5)
            sage: w._poly_rep()
            x
            sage: P = W(5)._poly_rep(); P
            ((a^2 + 4*a)*5^19 + (a^2 + a)*5^20)*x + 5 + (2*a^2 + 2*a + 3)*5^19 + (a^2 + 4*a)*5^20

        The coefficients of P are floating point p-adics::

            sage: ring = P.parent().base_ring()
            sage: ring
            5-adic Unramified Extension Ring in a defined by x^3 + 3*x + 3
            sage: ring._prec_type()
            'floating-point'
        """
        return self.unit.parent()(ccoefficients(self.unit, self.ordp, self.prime_pow.ram_prec_cap, self.prime_pow))
