
#*****************************************************************************
#   Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.categories.morphism import Morphism
from constructor import EllipticCurve
from sage.categories.homset import Hom

class WeierstrassIsomorphism(Morphism):

    def __init__(self, E, F):
        """
        Given two Weierstrass models $E$ and $F$ of the same elliptic curve,
        construct an isomorphism from $E$ to $F$.

        Alternatively, the second argument may be a vector of the form [u,r,s,t]
        which specifies a transformation

            $$ (x,y) \mapsto (x',y') = (u^2*x+r , u^3*y + s*u^2*x' + t) $$

        It is usually easier to use the \code{isomorphism_to} method of elliptic curves.
        """
        from ell_generic import is_EllipticCurve
        if is_EllipticCurve(E):
            K = F.base_ring()
            a1, a2, a3, a4, a6 = E.ainvs()
            b1, b2, b3, b4, b5 = F.ainvs()
            try:
                D = (F.discriminant() / E.discriminant())
                u = D.nth_root(12)
            except ValueError:
                raise ValueError, "Elliptic curves not isomorphic."
            if u.parent() is not K or K.is_exact() and u**12 != D:
                raise ValueError, "Elliptic curves not isomorphic."
            s = (a1*u - b1)/2
            r = (a2*u**2 + a1*s*u - s**2 - b2)/3
            t = (a3*u**3 - a1*r*u + 2*r*s - b3)/2
            urst = u, r, s, t
            if F != E.change_ring(K).change_weierstrass_model(urst):
                raise ValueError, "Elliptic curves not isomorphic."
        else:
            urst = F
            F = E.change_weierstrass_model(urst)
        Morphism.__init__(self, Hom(E(0).parent(), F(0).parent()))
        self._urst = urst
        self._domain_curve = E
        self._codomain_curve = F

    def _call_(self, P):
        if P[2] == 0:
            return self._codomain_curve(0)
        else:
            x, y, z = P
            u, r, s, t = self._urst
            new_x = u**2*x+r
            new_y = u**3*y + s*u**2*x + t
            return self._codomain_curve((new_x, new_y, 1))

    def __invert__(self):
        """
        EXAMPLES:
            sage: E = EllipticCurve('5077')
            sage: F = E.change_weierstrass_model([2,3,4,5]); F
            Elliptic Curve defined by y^2 - 8*x*y + 22*y = x^3 - 25*x^2 + 3*x + 588 over Rational Field
            sage: w = E.isomorphism_to(F)
            sage: P = E(-2,3,1)
            sage: w(P)
            (-5 : -3 : 1)
            sage: ~w
            Generic morphism:
            From: Abelian group of points on Elliptic Curve defined by y^2 - 8*x*y + 22*y = x^3 - 25*x^2 + 3*x + 588 over Rational Field
            To:   Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
            sage: Q = w(P); Q
            (-5 : -3 : 1)
            sage: (~w)(Q)
            (-2 : 3 : 1)
        """
        return WeierstrassIsomorphism(self._codomain_curve, self._domain_curve)

