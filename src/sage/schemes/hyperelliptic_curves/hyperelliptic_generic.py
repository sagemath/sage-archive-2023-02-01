"""
Hyperelliptic curves over a general ring

EXAMPLE::

    sage: P.<x> = GF(5)[]
    sage: f = x^5 - 3*x^4 - 2*x^3 + 6*x^2 + 3*x - 1
    sage: C = HyperellipticCurve(f); C
    Hyperelliptic Curve over Finite Field of size 5 defined by y^2 = x^5 + 2*x^4 + 3*x^3 + x^2 + 3*x + 4

EXAMPLE::

    sage: P.<x> = QQ[]
    sage: f = 4*x^5 - 30*x^3 + 45*x - 22
    sage: C = HyperellipticCurve(f); C
    Hyperelliptic Curve over Rational Field defined by y^2 = 4*x^5 - 30*x^3 + 45*x - 22
    sage: C.genus()
    2

    sage: D = C.affine_patch(0)
    sage: D.defining_polynomials()[0].parent()
    Multivariate Polynomial Ring in x0, x1 over Rational Field
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import PolynomialRing
from sage.misc.all import latex

import sage.schemes.plane_curves.projective_curve as plane_curve
from sage.schemes.generic.all import ProjectiveSpace

def is_HyperellipticCurve(C):
    """
    EXAMPLES::

        sage: R.<x> = QQ[]; C = HyperellipticCurve(x^3 + x - 1); C
        Hyperelliptic Curve over Rational Field defined by y^2 = x^3 + x - 1
        sage: sage.schemes.hyperelliptic_curves.hyperelliptic_generic.is_HyperellipticCurve(C)
        True
    """
    return isinstance(C,HyperellipticCurve_generic)

class HyperellipticCurve_generic(plane_curve.ProjectiveCurve_generic):
    def __init__(self, PP, f, h=None, names=None, genus=None):
        x, y, z = PP.gens()
        df = f.degree()
        F1 = sum([ f[i]*x**i*z**(df-i) for i in range(df+1) ])
        if h is None:
            F = y**2*z**(df-2) - F1
        else:
            dh = h.degree()
            deg = max(df,dh+1)
            F0 = sum([ h[i]*x**i*z**(dh-i) for i in range(dh+1) ])
            F = y**2*z**(deg-2) + F0*y*z**(deg-dh-1) - F1*z**(deg-df)
        plane_curve.ProjectiveCurve_generic.__init__(self,PP,F)
        R = PP.base_ring()
        if names == None:
            names = ["x","y"]
        elif isinstance(names,str):
            names = names.split(",")
        self._names = names
        P1 = PolynomialRing(R,name=names[0])
        P2 = PolynomialRing(P1,name=names[1])
        self._PP = PP
        self._printing_ring = P2
        self._hyperelliptic_polynomials = (f,h)
        self._genus = genus

    def change_ring(self, R):
        from constructor import HyperellipticCurve
        f, h = self._hyperelliptic_polynomials
        y = self._printing_ring.gen()
        x = self._printing_ring.base_ring().gen()
        return HyperellipticCurve(f.change_ring(R), h, "%s,%s"%(x,y))

    def _repr_(self):
        """
        String representation of hyperelliptic curves.

        EXAMPLE::

            sage: P.<x> = QQ[]
            sage: f = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C = HyperellipticCurve(f); C
            Hyperelliptic Curve over Rational Field defined by y^2 = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C = HyperellipticCurve(f,names='u,v'); C
            Hyperelliptic Curve over Rational Field defined by v^2 = 4*u^5 - 30*u^3 + 45*u - 22
        """

        f, h = self._hyperelliptic_polynomials
        R = self.base_ring()
        y = self._printing_ring.gen()
        x = self._printing_ring.base_ring().gen()
        if h == 0:
            return "Hyperelliptic Curve over %s defined by %s = %s"%(R, y**2, f(x))
        else:
            return "Hyperelliptic Curve over %s defined by %s + %s = %s"%(R, y**2, h(x)*y, f(x))

    def __cmp__(self, other):
        if not isinstance(other, HyperellipticCurve_generic):
            return -1
        return cmp(self._hyperelliptic_polynomials, other._hyperelliptic_polynomials)

    def hyperelliptic_polynomials(self, K=None, var='x'):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]; C = HyperellipticCurve(x^3 + x - 1, x^3/5); C
            Hyperelliptic Curve over Rational Field defined by y^2 + 1/5*x^3*y = x^3 + x - 1
            sage: C.hyperelliptic_polynomials()
            (x^3 + x - 1, 1/5*x^3)
        """
        if K is None:
            return self._hyperelliptic_polynomials
        else:
            f, h = self._hyperelliptic_polynomials
            P = PolynomialRing(K, var)
            return (P(f),P(h))

    def lift_x(self, x, all=False):
        f, h = self._hyperelliptic_polynomials
        x += self.base_ring()(0)
        one = x.parent()(1)
        if h.is_zero():
            y2 = f(x)
            if y2.is_square():
                if all:
                    return [self.point([x, y, one], check=False) for y in y2.sqrt(all=True)]
                else:
                    return self.point([x, y2.sqrt(), one], check=False)
        else:
            b = h(x)
            D = b*b + 4*f(x)
            if D.is_square():
                if all:
                    return [self.point([x, (-b+d)/2, one], check=False) for d in D.sqrt(all=True)]
                else:
                    return self.point([x, (-b+D.sqrt())/2, one], check=False)
        if all:
            return []
        else:
            raise ValueError, "No point with x-coordinate %s on %s"%(x, self)


    def genus(self):
        return self._genus

    def jacobian(self):
        import jacobian_generic
        return jacobian_generic.HyperellipticJacobian_generic(self)

    def odd_degree_model(self):
        r"""
        Return an odd degree model of self, or raise ValueError if one does not exist over the field of definition.

        EXAMPLES::

            sage: x = QQ['x'].gen()
            sage: H = HyperellipticCurve((x^2 + 2)*(x^2 + 3)*(x^2 + 5)); H
            Hyperelliptic Curve over Rational Field defined by y^2 = x^6 + 10*x^4 + 31*x^2 + 30
            sage: H.odd_degree_model()
            Traceback (most recent call last):
            ...
            ValueError: No odd degree model exists over field of definition

            sage: K2 = QuadraticField(-2, 'a')
            sage: Hp2 = H.change_ring(K2).odd_degree_model(); Hp2
            Hyperelliptic Curve over Number Field in a with defining polynomial x^2 + 2 defined by y^2 = 6*a*x^5 - 29*x^4 - 20*x^2 + 6*a*x + 1

            sage: K3 = QuadraticField(-3, 'b')
            sage: Hp3 = H.change_ring(QuadraticField(-3, 'b')).odd_degree_model(); Hp3
            Hyperelliptic Curve over Number Field in b with defining polynomial x^2 + 3 defined by y^2 = -4*b*x^5 - 14*x^4 - 20*b*x^3 - 35*x^2 + 6*b*x + 1

            Of course, Hp2 and Hp3 are isomorphic over the composite
            extension.  One consequence of this is that odd degree models
            reduced over "different" fields should have the same number of
            points on their reductions.  43 and 67 split completely in the
            compositum, so when we reduce we find:

            sage: P2 = K2.factor(43)[0][0]
            sage: P3 = K3.factor(43)[0][0]
            sage: Hp2.change_ring(K2.residue_field(P2)).frobenius_polynomial()
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849
            sage: Hp3.change_ring(K3.residue_field(P3)).frobenius_polynomial()
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849
            sage: H.change_ring(GF(43)).odd_degree_model().frobenius_polynomial()
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849

            sage: P2 = K2.factor(67)[0][0]
            sage: P3 = K3.factor(67)[0][0]
            sage: Hp2.change_ring(K2.residue_field(P2)).frobenius_polynomial()
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489
            sage: Hp3.change_ring(K3.residue_field(P3)).frobenius_polynomial()
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489
            sage: H.change_ring(GF(67)).odd_degree_model().frobenius_polynomial()
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489

        TESTS::
            sage: HyperellipticCurve(x^5 + 1, 1).odd_degree_model()
            Traceback (most recent call last):
            ...
            NotImplementedError: odd_degree_model only implemented for curves in Weierstrass form

            sage: HyperellipticCurve(x^5 + 1, names="U, V").odd_degree_model()
            Hyperelliptic Curve over Rational Field defined by V^2 = U^5 + 1
        """
        f, h = self._hyperelliptic_polynomials
        if h:
            raise NotImplementedError, "odd_degree_model only implemented for curves in Weierstrass form"
        if f.degree() % 2:
            # already odd, so just yield self
            return self

        rts = f.roots(multiplicities=False)
        if not rts:
            raise ValueError, "No odd degree model exists over field of definition"
        rt = rts[0]
        x = f.parent().gen()
        fnew =  f((x*rt + 1)/x).numerator() # move rt to "infinity"

        from constructor import HyperellipticCurve
        return HyperellipticCurve(fnew, 0, names=self._names, PP=self._PP)

    def has_odd_degree_model(self):
        r"""
        Return True if an odd degree model of self exists over the field of definition; False otherwise.

        Use ``odd_degree_model`` to calculate an odd degree model.

        EXAMPLES::
            sage: x = QQ['x'].0
            sage: HyperellipticCurve(x^5 + x).has_odd_degree_model()
            True
            sage: HyperellipticCurve(x^6 + x).has_odd_degree_model()
            True
            sage: HyperellipticCurve(x^6 + x + 1).has_odd_degree_model()
            False
        """
        try:
            return bool(self.odd_degree_model())
        except ValueError:
            return False

    def _magma_init_(self, magma):
        """
        Internal function. Returns a string to initialize this elliptic
        curve in the Magma subsystem.

        EXAMPLES::

            sage: R.<x> = QQ[]; C = HyperellipticCurve(x^3 + x - 1, x); C
            Hyperelliptic Curve over Rational Field defined by y^2 + x*y = x^3 + x - 1
            sage: magma(C)             # optional - magma
            Hyperelliptic Curve defined by y^2 + x*y = x^3 + x - 1 over Rational Field
            sage: R.<x> = GF(9,'a')[]; C = HyperellipticCurve(x^3 + x - 1, x^10); C
            Hyperelliptic Curve over Finite Field in a of size 3^2 defined by y^2 + x^10*y = x^3 + x + 2
            sage: D = magma(C); D             # optional - magma
            Hyperelliptic Curve defined by y^2 + (x^10)*y = x^3 + x + 2 over GF(3^2)
            sage: D.sage()                    # optional - magma
            Hyperelliptic Curve over Finite Field in a of size 3^2 defined by y^2 + x^10*y = x^3 + x + 2
        """
        f, h = self._hyperelliptic_polynomials
        return 'HyperellipticCurve(%s, %s)'%(f._magma_init_(magma), h._magma_init_(magma))


    def monsky_washnitzer_gens(self):
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        return S.gens()
