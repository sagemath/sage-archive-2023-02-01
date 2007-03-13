r"""
p-Adic Extension Rings

DESCRIPTION:

Extensions of $\Z_p$ are implemented as polynomials over $\Z_p$.
Given a polynomial to generate the extension, we use the Round4
Algorithm to find an unramified and totally ramified polynomial
which together generate the required extension.

"""


from sage.rings.padics.padic_ring_generic import pAdicRingGeneric
from sage.rings.polynomial_quotient_ring import *


class pAdicRingExtension(pAdicRingGeneric, PolynomialQuotientRing_domain):
    r"""
    Base class for extensions of Z_p.

    As of right now, implemented in the simplest way possible:
    we just take a polynomial and a ring, and return the polynomial
    quotient ring, after checking that the defining polynomial is
    irreducible.

    TODO:
    We should use the Round4 algorithm to break the given
    polynomial into a totally ramified and unramified part, and
    create two subextensions, so that the whole extension is
    given as a tower.
    """

    def __init__(self, poly):
        R = poly.base_ring()
        pAdicRingGeneric.__init__(self, R.prime(), R.precision_cap(), R.print_mode())
        PolynomialQuotientRing_domain.__init__(self, poly.parent(), poly)

    def __call__(self, x):
        return PolynomialQuotientRing_domain.__call__(self, x)
#        r"""
#        This is currently unimplemented, as understanding how to
#        embed one ring in another when appropriate is nontrivial.
#        """
#        raise NotImplementedError, "not smart enough to embed rings"

    def _repr_(self):
        return "Extension of %s in %s defined by %s"%(
            self.base_ring(), self.variable_name(), self.modulus())

    def fraction_field(self):
        r"""
        Returns the fraction field of this extension, which is just
        the extension of base.fraction_field() determined by the
        same polynomial.
        """

        return self.base_ring().fraction_field().extension(self._poly)


