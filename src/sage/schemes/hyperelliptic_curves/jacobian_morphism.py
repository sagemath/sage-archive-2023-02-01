"""
Jacobian ``morphism'' as a class in the Picard group
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.generic.morphism import SchemeMorphism

from sage.rings.all import PolynomialRing, ZZ

def cantor_reduction_simple(a1,b1,f,genus):
    # Divisor reduction.
    a2 = (f - b1**2)//a1
    a2 *= 1/a2.leading_coefficient()
    b2 = -b1.mod(a2);
    if a2.degree() == a1.degree():
        assert a2.degree() == genus+1
        print "Returning ambiguous form of degree genus+1."
        return (a2, b2)
    elif a2.degree() > genus:
        return cantor_reduction_simple(a2,b2,f)
    return (a2, b2)

def cantor_reduction(a,b,f,h,genus):
    assert a.degree() < 2*genus+1
    assert b.degree() < a.degree()
    k = f - h*b - b**2
    if 2*a.degree() == k.degree():
        # must adjust b to include the point at infinity
        g1 = a.degree()
        x = a.parent().gen()
        r = (x**2 + h[g1]*x - f[2*g1]).roots()[0][0]
        b = b + r*(x**g1 - (x**g1).mod(a))
        k = f - h*b - b**2
    assert k.mod(a) == 0
    a = k//a
    a //= a.leading_coefficient()
    b = -(b+h).mod(a)
    if a.degree() > genus:
        return cantor_reduction(a,b,f,h,genus)
    return (a, b)

def cantor_composition_simple(D1,D2,f,genus):
    a1, b1 = D1
    a2, b2 = D2
    if a1 == a2 and b1 == b2:
        # Duplication law:
        d, h1, h3 = a1.xgcd(2*b1)
        a = (a1//d)**2
        b = (b1 + h3*((f - b1**2)//d)).mod(a)
    else:
        d0, _, h2 = a1.xgcd(a2)
        if d0 == 1:
            a = a1*a2
            b = (b2 + h2*a2*(b1-b2)).mod(a)
        else:
            d, l, h3 = d0.xgcd(b1 + b2)
            a = (a1*a2)//(d**2)
            b = (b2 + l*h2*(b1-b2)*(a2//d) + h3*((f - b2**2)//d)).mod(a)
    if a.degree() > genus:
        return cantor_reduction_simple(a,b,f,genus)
    return (a,b)

def cantor_composition(D1,D2,f,h,genus):
    a1, b1 = D1
    a2, b2 = D2
    if a1 == a2 and b1 == b2:
        # Duplication law:
        d, h1, h3 = a1.xgcd(2*b1 + h)
        # NOTE THAT d is not normalised, but this gives a crash:
        # d *= 1/d.leading_coefficient()
        # print "d =", d
        a = (a1//d)**2;
        b = (b1 + h3*((f-h*b1-b1**2)//d)).mod(a)
    else:
        d0, _, h2 = a1.xgcd(a2)
        if d0 == 1:
            a = a1*a2;
            b = (b2 + h2*a2*(b1-b2)).mod(a)
        else:
            e0 = b1+b2+h
            if e0 == 0:
                a = (a1*a2)//(d0**2);
                b = (b2 + h2*(b1-b2)*(a2//d0)).mod(a)
            else:
                d, l, h3 = d0.xgcd(e0)
                a = (a1*a2)//(d**2);
                b = (b2 + l*h2*(b1-b2)*(a2//d) + h3*((f-h*b2-b2**2)//d)).mod(a)
    a *= 1/a.leading_coefficient()
    if a.degree() > genus:
        return cantor_reduction(a,b,f,h,genus)
    return (a,b)

class JacobianMorphism_divisor_class(SchemeMorphism):
    r"""
    An element of a $J(K) = \Pic^0_K(C)$.
    """
    def __init__(self, parent, polys, reduce=True, check=False):
        SchemeMorphism.__init__(self, parent)
        if polys == 0:
            P = PolynomialRing(self.parent(), 'x')
            self.__polys = (P(1),P(0))
        if check:
            C = parent.curve()
            K = parent.value_ring()
            f, h = C.hyperelliptic_polynomials(K)
            a, b = polys
            if not (b**2 + h*b - f)%a == 0:
                raise ValueError, \
                      "Argument polys (= %s) must be reduced divisor on curve %s."%(
                    polys, C)
        self.__polys = polys

    def __repr__(self):
        a, b = self.__polys
        if a == 1:
            return "(1)"
        P = self.parent()._printing_ring
        y = P.gen()
        x = P.base_ring().gen()
        return "(%s, %s)"%(a(x), y - b(x))

    def list(self):
        return self.__polys

    def __add__(self,other):
        X = self.parent()
        C = X.curve()
        K = X.value_ring()
        f, h = C.hyperelliptic_polynomials(K)
        if h == 0:
            D = cantor_composition_simple(self.__polys,other.__polys,f,C.genus())
        else:
            D = cantor_composition(self.__polys,other.__polys,f,h,C.genus())
        return JacobianMorphism_divisor_class(X, D, reduce=False, check=False)

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        if self.is_zero():
            return self
        polys = self.__polys
        X = self.parent()
        K = X.value_ring()
        f, h = X.curve().hyperelliptic_polynomials(K)
        if h == 0:
            D = (polys[0],-polys[1])
        else:
            D = (polys[0],-polys[1]-h.mod(polys[0]))
        return JacobianMorphism_divisor_class(X, D, reduce=False, check=False)

    def __mul__(self, n):
        try:
            n = ZZ(n)
        except TypeError:
            raise TypeError, "Argument n (= %s) must be an integer."
        X = self.parent()
        if n < 0:
            return self * (-n)
        elif n == 0:
            P = PolynomialRing(X.value_ring(), 'x')
            D = (P(1),P(0))
            return X(0)
        elif n == 1:
            return self
        m = n//2
        return self * m + self * (n-m)

    def __lmul__(self, n):
        return self.__mul__(n)

    def __rmul__(self, n):
        return self.__mul__(n)

    def __nonzero__(self):
        return self.__polys[0] != 1
