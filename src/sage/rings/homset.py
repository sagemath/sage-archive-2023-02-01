"""
Space of homomorphisms between two rings.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.all import HomsetWithBase, Rings
from sage.structure.parent_base import ParentWithBase

import morphism
import quotient_ring

RINGS = Rings()


def is_RingHomset(H):
    return isinstance(H, RingHomset_generic)

def RingHomset(R, S):
    if quotient_ring.is_QuotientRing(R):
        return RingHomset_quo_ring(R, S)
    return RingHomset_generic(R, S)


class RingHomset_generic(HomsetWithBase):
    def __init__(self, R, S):
        HomsetWithBase.__init__(self, R, S, RINGS)

    def _repr_(self):
        return "Set of Homomorphisms from %s to %s"%(self.domain(), self.codomain())

    def has_coerce_map_from(self, x):
        """
        The default for coercion maps between ring homomorphism
        spaces is very restrictive (until more implementation work
        is done).
        """
        return (x.domain() == self.domain() and x.codomain() == self.codomain())

    def _coerce_impl(self, x):
        if not isinstance(x, morphism.RingHomomorphism):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            if isinstance(x, morphism.RingHomomorphism_im_gens):
                return morphism.RingHomomorphism_im_gens(self, x.im_gens())
            elif isinstance(x, morphism.RingHomomorphism_cover):
                return morphism.RingHomomorphism_cover(self)
        raise TypeError

    def __call__(self, im_gens, check=True):
        """
        EXAMPLES:
            sage: H = Hom(ZZ, QQ)
            sage: phi = H([])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism

        TESTS:
            sage: H = Hom(ZZ, QQ)
            sage: H == loads(dumps(H))
            True
        """
        if isinstance(im_gens, (morphism.RingHomomorphism_im_gens,  morphism.RingHomomorphism_cover) ):
            return self._coerce_impl(im_gens)
        try:
            return morphism.RingHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError), err:
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError, "images do not define a valid homomorphism"


    def natural_map(self):
        return morphism.RingHomomorphism_coercion(self)


class RingHomset_quo_ring(RingHomset_generic):
    """
    Space of ring homomorphism where the domain is a (formal) quotient ring.

    EXAMPLES:
        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quotient(x^2 + y^2)
        sage: phi = S.hom([b,a]); phi
        Ring endomorphism of Quotient of Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
          Defn: a |--> b
                b |--> a
        sage: phi(a)
        b
        sage: phi(b)
        a

    TESTS:
    We test pickling of a homset from a quotient.
        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quotient(x^2 + y^2)
        sage: H = S.Hom(R)
        sage: H == loads(dumps(H))
        True

    We test pickling of actual homomorphisms in a quotient:
        sage: phi = S.hom([b,a])
        sage: phi == loads(dumps(phi))
        True
    """
    def __call__(self, im_gens, check=True):
        if isinstance(im_gens, morphism.RingHomomorphism_from_quotient):
            return morphism.RingHomomorphism_from_quotient(self, im_gens._phi())
        try:
            pi = self.domain().cover()
            phi = pi.domain().hom(im_gens, check=check)
            return morphism.RingHomomorphism_from_quotient(self, phi)
        except (NotImplementedError, ValueError), err:
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError, "images do not define a valid homomorphism"

    def _coerce_impl(self, x):
        if not isinstance(x, morphism.RingHomomorphism_from_quotient):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return morphism.RingHomomorphism_from_quotient(self, x._phi())
        raise TypeError
