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

    def __call__(self, im_gens, check=True):
        """
        EXAMPLES:
            sage: H = Hom(ZZ, QQ)
            sage: phi = H([])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism
        """
        try:
            return morphism.RingHomomorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError), err:
            raise TypeError, "%s\nimages do not define a valid homomorphism"%err

    def natural_map(self):
        return morphism.RingHomomorphism_coercion(self)


class RingHomset_quo_ring(RingHomset_generic):
    """
    Space of ring homomorphism where the domain is a (formal) quotient ring.

    EXAMPLES:
        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quotient(x^2 + y^2)
        sage: phi = S.hom([b,a]); phi
        Ring endomorphism of Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y^2 + x^2)
          Defn: a |--> b
                b |--> a
        sage: phi(a)
        b
        sage: phi(b)
        a
    """
    def __call__(self, im_gens, check=True):
        try:
            pi = self.domain().cover()
            phi = pi.domain().hom(im_gens, check=check)
            return morphism.RingHomomorphism_from_quotient(self, phi)
        except (NotImplementedError, ValueError), err:
            raise TypeError, "images do not define a valid homomorphism"


