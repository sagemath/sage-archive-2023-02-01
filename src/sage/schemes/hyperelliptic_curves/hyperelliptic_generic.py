"""
Hyperelliptic curves over a general ring

EXAMPLE:
    sage: P, x = PolynomialRing(GF(5),"x").objgen()
    sage: f = x^5 - 3*x^4 - 2*x^3 + 6*x^2 + 3*x - 1
    sage: C = HyperellipticCurve(f); C
    Hyperelliptic Curve over Finite Field of size 5 defined by y^2 = x^5 + 2*x^4 + 3*x^3 + x^2 + 3*x + 4

EXAMPLE:
    sage: P, x = PolynomialRing(QQ,"x").objgen()
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
    EXAMPLES:
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
        P1 = PolynomialRing(R,name=names[0])
        P2 = PolynomialRing(P1,name=names[1])
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
        String representation hyperelliptic curves.

        EXAMPLE:
            sage: P, x = PolynomialRing(QQ,"x").objgen()
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
        EXAMPLES:
            sage: R.<x> = QQ[]; C = HyperellipticCurve(x^3 + x - 1, x^3/5); C
            Hyperelliptic Curve over Rational Field defined by y^2 + 1/5*x^3*y = x^3 + x - 1
            sage: C.hyperelliptic_polynomials()
            (x^3 + x - 1, 1/5*x^3)
        """
        if K == None:
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

    def _magma_init_(self, magma):
        """
        Internal function.  Returns a string to initialize this
        elliptic curve in the Magma subsystem.

        EXAMPLES:
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
