from __future__ import print_function, absolute_import, division

from sage.rings.derivation import RingDerivation
from sage.structure.category_object import normalize_names
from sage.categories.morphism import Morphism, IdentityMorphism
import sage.rings.ring as ring
from sage.rings.finite_rings.finite_field_base import is_FiniteField



def OrePolynomialRing(base_ring, base_ring_automorphism = None, base_ring_derivation = None, names = None, sparse = False):
    if base_ring not in sage.categories.rings.Rings().Commutative():
        raise TypeError("base_ring must be a commutative ring")
    if base_ring_automorphism is None:
        base_ring_automorphism = IdentityMorphism(base_ring)
    else:
        if (not isinstance(base_ring_automorphism, Morphism)
            or base_ring_automorphism.domain() != base_ring
            or base_ring_automorphism.codomain() != base_ring):
            raise TypeError("base_ring_automorphism must be a ring automorphism of base_ring (=%s)" % base_ring)
    if base_ring_derivation is None:
        if not isinstance(base_ring, (PolynomialRing_general, MpolynomialRing_generic)):
            raise NotImplementedError()
        else:
            H = Homset(base_ring, base_ring)
            base_ring_derivation = RingDerivation_polynomial(H, base_ring_automorphisme, [0 for _ in range(base_ring.ngens())])
    else:
        if (not isinstance(base_ring_derivation, RingDerivation)
            or base_ring_derivation.domain() != base_ring
            or base_ring_derivation.codomain() != base_ring):
            raise TypeError("base_ring_derivation must be a derivation of base_ring (=%s)" % base_ring)
        elif base_ring_derivation._theta != base_ring_automorphism:
            raise TypeError("base_ring_derivation must be a base_ring_automorphism-derivation")
    if sparse:
        raise NotImplementedError("Sparse Ore polynomial rings are not implemented")
    if names is None:
        raise TypeError("you must specify the name of the variable")
    try:
        names = normalize_names(1, names)[0]
    except IndexError:
        raise NotImplementedError("multivariate Ore polynomials rings not supported")
    from sage.rings.polynomial.Ore_polynomial_ring import OrePolynomialRing_general
    return OrePolynomialRing_general(base_ring, base_ring_automorphism, base_ring_derivation, names, sparse)
    


def OrePolynomialRing1(base_ring, base_ring_automorphism = None, base_ring_derivation = None, names = None, sparse = False):
    if names is None:
        raise TypeError("you must specify the name of the variable")
    try:
        names = normalize_names(1, names)[0]
    except IndexError:
        raise NotImplementedError("multivariate Ore polynomials rings not supported")
    from sage.rings.polynomial.Ore_polynomial_ring import OrePolynomialRing_general
    return OrePolynomialRing_general(base_ring, base_ring_automorphism, base_ring_derivation, names, sparse)
