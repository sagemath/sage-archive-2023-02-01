"""
Rational point sets on a Jacobian

EXAMPLES::

    sage: x = QQ['x'].0
    sage: f = x^5 + x + 1
    sage: C = HyperellipticCurve(f); C
    Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1
    sage: C(QQ)
    Set of Rational Points of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1
    sage: P = C([0,1,1])
    sage: J = C.jacobian(); J
    Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1
    sage: Q = J(QQ)(P); Q
    (x, y - 1)
    sage: Q + Q
    (x^2, y - 1/2*x - 1)
    sage: Q*3
    (x^2 - 1/64*x + 1/8, y + 255/512*x + 65/64)
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.schemes.generic.spec as spec
from sage.rings.all import is_Polynomial, PolynomialRing, Integer, ZZ
from sage.schemes.generic.homset import SchemeHomset_generic
from sage.schemes.generic.morphism import is_SchemeMorphism
#from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic
from hyperelliptic_generic import is_HyperellipticCurve
from jacobian_morphism import JacobianMorphism_divisor_class_field

class JacobianHomset_divisor_classes(SchemeHomset_generic):
    def __init__(self, X, S):
        R = X.base_ring()
        SchemeHomset_generic.__init__(self, spec.Spec(S, R), X)
        P2 = X.curve()._printing_ring
        if S != R:
            y = str(P2.gen())
            x = str(P2.base_ring().gen())
            P1 = PolynomialRing(S,name=x)
            P2 = PolynomialRing(P1,name=y)
        self._printing_ring = P2

    def __call__(self, P):
        r"""
        Returns a rational point P in the abstract Homset J(K), given:

        0. A point P in J = Jac(C), returning P; 1. A point P on the curve
        C such that J = Jac(C), where C is an odd degree model, returning
        [P - oo]; 2. A pair of points (P, Q) on the curve C such that J =
        Jac(C), returning [P-Q]; 2. A list of polynomials (a,b) such that
        `b^2 + h*b - f = 0 mod a`, returning [(a(x),y-b(x))].

        EXAMPLES::

            sage: P.<x> = PolynomialRing(QQ)
            sage: f = x^5 - x + 1; h = x
            sage: C = HyperellipticCurve(f,h,'u,v')
            sage: P = C(0,1,1)
            sage: J = C.jacobian()
            sage: Q = J(QQ)(P)
            sage: for i in range(6): i*Q
            (1)
            (u, v - 1)
            (u^2, v + u - 1)
            (u^2, v + 1)
            (u, v + 1)
            (1)
        """
        if isinstance(P,(int,long,Integer)) and P == 0:
            R = PolynomialRing(self.value_ring(), 'x')
            return JacobianMorphism_divisor_class_field(self, (R(1),R(0)))
        elif isinstance(P,(list,tuple)):
            if len(P) == 1 and P[0] == 0:
                R = PolynomialRing(self.value_ring(), 'x')
                return JacobianMorphism_divisor_class_field(self, (R(1),R(0)))
            elif len(P) == 2:
	        P1 = P[0]; P2 = P[1]
		if is_Polynomial(P1) and is_Polynomial(P2):
                    return JacobianMorphism_divisor_class_field(self, tuple(P))
                if is_SchemeMorphism(P1) and is_SchemeMorphism(P2):
                    return self(P1) - self(P2)
            raise TypeError, "Argument P (= %s) must have length 2."%P
        elif isinstance(P,JacobianMorphism_divisor_class_field) and self == P.parent():
            return P
        elif is_SchemeMorphism(P):
            x0 = P[0]; y0 = P[1]
            R, x = PolynomialRing(self.value_ring(), 'x').objgen()
            return self((x-x0,R(y0)))
        raise TypeError, "Argument P (= %s) does not determine a divisor class"%P

    def _cmp_(self,other):
        if self.curve() == other.curve():
            return 0
        else:
            return -1

    def _point_morphism_class(self, *args, **kwds):
        return JacobianMorphism_divisor_class_field(*args, **kwds)

    def curve(self):
        return self.codomain().curve()

    def value_ring(self):
        """
        Returns S for a homset X(T) where T = Spec(S).
        """
        T = self.domain()
        if spec.is_Spec(T):
            return T.coordinate_ring()
        else:
            raise TypeError, "Domain of argument must be of the form Spec(S)."

    def base_extend(self, R):
        if R != ZZ:
            raise NotImplementedError, "Jacobian point sets viewed as modules over rings other than ZZ not implemented"
        return self
