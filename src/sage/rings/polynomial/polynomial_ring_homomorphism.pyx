"""
Ring homomorphisms from a polynomial ring to another ring

This module currently implements the canonical ring homomorphism from
`A[x]` to `B[x]` induced by a ring homomorphism from `A` to `B`.

.. TODO::

    Implement homomorphisms from `A[x]` to an arbitrary ring `R`,
    given by a ring homomorphism from `A` to `R` and the image of `x`
    in `R`.

AUTHORS:

- Peter Bruin (March 2014): initial version

"""

from sage.rings.morphism cimport RingHomomorphism_from_base
from sage.structure.element cimport Element

cdef class PolynomialRingHomomorphism_from_base(RingHomomorphism_from_base):
    """
    The canonical ring homomorphism from `R[x]` to `S[x]` induced by a
    ring homomorphism from `R` to `S`.

    EXAMPLE::

        sage: QQ['x'].coerce_map_from(ZZ['x'])
        Ring morphism:
          From: Univariate Polynomial Ring in x over Integer Ring
          To:   Univariate Polynomial Ring in x over Rational Field
          Defn: Induced from base ring by
                Natural morphism:
                  From: Integer Ring
                  To:   Rational Field

    """
    cpdef Element _call_(self, x):
        """
        Evaluate the homomorphism ``self`` at ``x``.

        EXAMPLE::

            sage: from sage.rings.polynomial.polynomial_ring_homomorphism import PolynomialRingHomomorphism_from_base
            sage: R.<x> = ZZ[]
            sage: S.<x> = QQ[]
            sage: f = ZZ.hom(QQ)
            sage: F = PolynomialRingHomomorphism_from_base(R.Hom(S), f)
            sage: F(2*x)
            2*x

        """
        P = self.codomain()
        f = self.underlying_map()
        return P([f(b) for b in x])

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        Evaluate ``self`` at ``x`` with additional (keyword) arguments.

        EXAMPLE::

            sage: from sage.rings.polynomial.polynomial_ring_homomorphism import PolynomialRingHomomorphism_from_base
            sage: R.<x> = ZZ[]
            sage: S.<x> = GF(5)[]
            sage: f = ZZ.hom(GF(5))
            sage: F = PolynomialRingHomomorphism_from_base(R.Hom(S), f)
            sage: F(2*x, check=True)
            2*x

        """
        P = self.codomain()
        f = self.underlying_map()
        return P([f(b) for b in x], *args, **kwds)
