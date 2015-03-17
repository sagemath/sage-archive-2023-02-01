r"""
AUTHORS:

- Ben Hutz (11-2012)
- Joao Alberto de Faria (10-2013)

"""
import sys
from sage.calculus.functions       import jacobian
from sage.categories.fields        import Fields
from sage.categories.number_fields import NumberFields
from sage.categories.homset        import Hom
from sage.functions.all            import sqrt
from sage.misc.cachefunc           import cached_method
from sage.misc.mrange              import mrange
from sage.rings.all                import Integer
from sage.rings.commutative_ring   import is_CommutativeRing
from sage.rings.finite_rings.constructor import GF
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.rings.fraction_field     import FractionField
from sage.rings.integer_ring       import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field     import QQ
from sage.rings.real_mpfr          import RealField
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme, AlgebraicScheme_subscheme_projective
from sage.schemes.product_projective.space import ProductProjectiveSpaces
from sage.symbolic.constants       import e
from copy                          import copy
_Fields = Fields()

def WehlerK3Surface(polys):
    r"""
    Defines a K3 Surface over $\mathbb{P}^2 \times \mathjbb{P}^2$ defined as
    the intersection of a bi-linear and bi-quadratic form.

    INPUT: 
        Bilinear and Biquadratic polynomials as a tuple or list

    OUTPUT: 
        WhelerK3 Surface

    EXAMPLES:::

        sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
        sage: L = x0*y0+x1*y1-x2*y2
        sage: Q = x0*x1*y1^2 + x2^2*y0*y2
        sage: WehlerK3Surface([L,Q])
        Closed subscheme of Product of projective spaces P^2 x P^2 over Rational
        Field defined by:
         x0*y0 + x1*y1 - x2*y2,
         x0*x1*y1^2 + x2^2*y0*y2
    """

    R = polys[0].parent().base_ring()
    if R in _Fields:
        if is_FiniteField(R):
            return WehlerK3Surface_finite_field(polys)
        else:
            return WehlerK3Surface_field(polys)
    elif is_CommutativeRing(R):
        return WehlerK3Surface_ring(polys)
    else:
        raise TypeError, "R ( = %s) must be a commutative ring"%R

def random_WehlerK3Surface(PP):
    r"""

    Produces a random K3 surface in $\mathbb{P}^2 \times \mathjbb{P}^2$ defined as the
    intersection of a bi-linear and bi-quadratic form.

    INPUT: Projective Space Cartesian Product

    OUTPUT: WehlerK3 Surface
    """
    def getmon():
        return ([x, y] for x in range(0,3) for y in range(0,3))

    CR = PP.coordinate_ring()
    BR = PP.base_ring()
    Q = 0
    for a in getmon():
        for b in getmon():
            Q += BR.random_element() * CR.gen(a[0]) * CR.gen(a[1]) * CR.gen(3+b[0]) * CR.gen(3+b[1])
    L = CR.gen(0) * CR.gen(3) + CR.gen(1) * CR.gen(4) + CR.gen(2) * CR.gen(5)
    return( WehlerK3Surface([L,Q]) )

class WehlerK3Surface_ring(AlgebraicScheme_subscheme_projective):
    r"""
    #TODO - citations

    A K3 surface in $\mathbb{P}^2 \times \mathjbb{P}^2$ defined as the
    intersection of a bi-linear and bi-quadratic form.

    EXAMPLES:::

        sage: R.<x,y,z,u,v,w> = PolynomialRing(QQ,6)
        sage: L = x*u-y*v
        sage: Q = x*y*v^2 + z^2*u*w
        sage: WehlerK3Surface([L,Q])
        Closed subscheme of Product of projective spaces P^2 x P^2 over Rational
        Field defined by:
          x*u - y*v,
          x*y*v^2 + z^2*u*w
    """
    def __init__(self, polys):
        R = polys[0].parent()
        vars = R.variable_names()
        A = ProductProjectiveSpaces([ZZ(2),ZZ(2)],R.base_ring(),vars)
        CR = A.coordinate_ring()
        """Check for following:
            Is the user calling in 2 polynomials from a list or tuple?
            Is there one biquadratic and one bilinear polynomial?
        """
        if not isinstance(polys, (list,tuple)):
            raise TypeError, "polys must be a list or tuple of polynomials"
        if len(polys) != 2:
            raise AttributeError, "There must be 2 polynomials"

        if ( all ((( e[0] + e[1] + e[2]) == 1 and ( e[3] + e[4] + e[5]) == 1) for e in polys[0].exponents() )):
            self.L = CR(polys[0])
        elif ( all ((( e[0] + e[1] + e[2]) == 1 and ( e[3] + e[4] + e[5]) == 1) for e in polys[1].exponents() )):
            self.L = CR(polys[1])
        else:
            raise AttributeError, "There must be one bilinear polynomial"

        if ( all ((( e[0] + e[1] + e[2]) == 2 and ( e[3] + e[4] + e[5]) == 2) for e in polys[0].exponents() )):
            self.Q = CR(polys[0])
        elif ( all ((( e[0] + e[1] + e[2]) == 2 and ( e[3] + e[4] + e[5]) == 2) for e in polys[1].exponents() )):
            self.Q = CR(polys[1])
        else:
            raise AttributeError, "There must be one biquadratic polynomial"
        AlgebraicScheme_subscheme.__init__(self, A, polys)
    def _morphism(self, *args, **kwds):
        return ProductProjectiveSpaces_morphism(*args, **kwds)

    def change_ring(self,R):
        r"""
        Changes the base ring on which the WhelerK3 Surface is defined

        INPUT:
            'R' - ring
        OUTPUT: 
            WhelerK3 Surface defined over input ring
        EXAMPLES:::
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],GF(3))
            sage: L = x0*y0+x1*y1-x2*y2
            sage: Q = x0*x1*y1^2 + x2^2*y0*y2
            W = WehlerK3Surface([L,Q])
            W.base_ring()  #Finite Field of size 3
            T = W.change_ring(GF(7))
            T.base_ring()
        """
        if R.is_ring() == False:
            return TypeError, "R must be a ring"

        LR = self.L.change_ring(R)
        LQ = self.Q.change_ring(R)
        return (WehlerK3Surface( [LR,LQ]))

    def _check_satisfies_equations(self, P):
        r"""
        Function checks to see if point lies on the Wehler K3 Surface

        INPUT:
            'P' - point in $\mathbb{P}^2 \times \mathbb{P}^2$
        OUTPUT:
            Error if the point is not on the surface
        EXAMPLES:::
            sage: P.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: P = X([0,0,1,1,0,0])
            sage:X._check_satisfies_equations(P)

            ::

            sage: P.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: P = X([0,1,1,1,0,0])
            sage:X._check_satisfies_equations(P)
            Traceback (click to the left of this block for traceback)
            ...
            TypeError: Fiber is degenerate
        """
        Point = list(P)
        if self.L(*Point) == 0 and self.Q(*Point) == 0:
            pass
        else:
            raise AttributeError, "Point must be on Surface"

    def _Lcoeff(self, component, i):
        r"""
        Returns the polynomials defined as:
            $L^x_i$ = the coefficients of $y_i$ in $L(x,y)$ ( Component = 0)
            $L^y_i$ = the coefficients of $x_i$ in $L(x,y)$ ( Component = 1)
            Definition and Notation borrowed from the following paper:
            G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.

        INPUT:
            'component' - Integer: 0 or 1
            'i' - Integer: 0, 1 or 2

        OUTPUT:
            Polynomial in terms of either y ( Component = 0) or x ( Component = 1)

        EXAMPLES::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X._Lcoeff(0,0)
             y0

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X._Lcoeff(1,0)
             x0
        """

        #Error Checks for Passed in Values
        if not component in [0,1]:
            raise ValueError, "Component can only be 1 or 0"
        if not i in [0,1,2]:
            raise ValueError, "Index must be 0, 1, or 2"
        R = self.ambient_space().coordinate_ring()
        return self.L.coefficient(R.gen(component*3 + i))

    def _Qcoeff(self, component, i, j):
        r"""
        Returns the polynomials defined as:
            $Q^x_{ij}$ = the coefficients of $y_{i}y_{j}$ in $Q(x,y)$ ( Component = 0)
            $Q^y_{ij}$ = the coefficients of $x_{i}x_{j}$ in $Q(x,y)$ ( Component = 1)
            Definition and notation borrowed from Call & Silverman

        INPUT:
            'component' - Integer: 0 or 1
            'i' - Integer: 0, 1 or 2
            'j' - Integer: 0, 1 or 2
        OUTPUT:
            Polynomial in terms of either y ( Component = 0) or x ( Component = 1)

        EXAMPLES::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X._Qcoeff(0,0,0)
             y0*y1 + y2^2

            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X._Qcoeff(1,1,0)
             x0^2

        """
        #Check Errors in Passed in Values
        if not component in [0,1]:
            raise ValueError, "Component can only be 1 or 0"

        if not (i in [0,1,2]) and  not (j in [0,1,2]):
            raise ValueError, "The two indexes must be either 0, 1, or 2"

        R = self.ambient_space().coordinate_ring()
        return self.Q.coefficient( R.gen( component * 3 + i) * R.gen( component * 3 + j))

    @cached_method
    def Gpoly( self, component, k):
        r"""
        Returns the G polynomials defined by $G^*_k = \left(L^*_j\right)^2Q^*_{ii}-L^*_iL^*_jQ^*_{ij}+\left(L^*_i\right)^2Q^*_{jj} $
        where {i,j,k} is some permutation of (0,1,2) and * is either x ( Component = 1) or y ( Component = 0)

        INPUT:
            'component' - Integer: 0 or 1
            'k' - Integer: 0, 1 or 2

        OUTPUT:
            Polynomial in terms of either y ( Component = 0) or x ( Component = 1)

        EXAMPLES::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.Gpoly(1,0)
             x0^2*x1^2 + x1^4 - x0*x1^2*x2 + x1^3*x2 + x1^2*x2^2 + x2^4
        """
        #Check Errors in passed in values
        if not component in [0,1]:
            raise ValueError, "Component can only be 1 or 0"

        if not k in [0,1,2]:
            raise ValueError, "Index must be either 0, 1, or 2"

        Indices = [ 0, 1, 2]
        Indices.remove( k)
        i = Indices[0]
        j = Indices[1]

        return (self._Lcoeff( component, j) ** 2) * ( self._Qcoeff( component, i, i) ) - ( self._Lcoeff( component, i) ) * \
            (self._Lcoeff( component, j) ) * (self._Qcoeff( component, i, j) ) + ( self._Lcoeff( component, i) ** 2) * \
            ( self._Qcoeff( component, j, j))

    @cached_method
    def Hpoly(self, component, i, j):
        r"""
        Returns the H polynomials defined by $H^*_{ij} = 2L^*_iL^*_jQ^*_{kk}-L^*_iL^*_kQ^*_{jk}--L^*_jL^*_kQ^*_{ik}+\left(L^*_k\right)^2Q^*_{ij}$
        where {i,j,k} is some permutation of (0,1,2) and * is either y ( Component = 0) or x ( Component = 1)

        INPUT:
            'component' - Integer: 0 or 1
            'i' - Integer: 0, 1 or 2
            'j' - Integer: 0, 1 or 2

        OUTPUT:
            Polynomial in terms of either y ( Component = 0) or x ( Component = 1)

        EXAMPLES::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.Hpoly(0,1,0)
             2*y0*y1^3 + 2*y0*y1*y2^2 - y1*y2^3

        """
        #Check Errors in Passed in Values
        if ( component in [0,1]) == False:
            raise ValueError, "Component can only be 1 or 0"

        if ( i in [0,1,2] == False) or ( j in [0,1,2] == False):
            raise ValueError, "The two indexes must be either 0, 1, or 2"

        Indices = [ 0, 1, 2]
        Indices.remove(i)
        Indices.remove(j)

        k = Indices[0]

        return 2*( self._Lcoeff( component, i)) * (self._Lcoeff( component, j)) * (self._Qcoeff( component, k, k)) - (self._Lcoeff( component, i)) * (self._Lcoeff( component, k)) * (self._Qcoeff( component, j, k)) - (self._Lcoeff( component, j)) * (self._Lcoeff( component, k)) * (self._Qcoeff( component, i, k)) + (self._Lcoeff( component, k) ** 2) * (self._Qcoeff( component, i, j))

    def Lxa(self,a):
        r"""
        Function will return a fiber defined by:
        $L^{x}_{a} = \{(a,y) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon L(a,y) = 0\}$
        Notation and definition from:
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.
        INPUT:
            'a' - Point in $\mathbb{P}^2$
        OUTPUT:
            A polynomial representing the fiber

        EXAMPLES::
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: T = PP(1,1,0,1,0,0);
            sage: X.Lxa(T[0])
              y0 + y1
        """
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[1]
        PSYC = PSY.coordinate_ring()
        p = ASC.hom([a[0],a[1],a[2]] +list(PSY.gens()),PSYC)
        return((p(self.L)))

    def Qxa(self,a):
        r"""
        Function will return a fiber defined by:
        $Q^{x}_{a} = \{(a,y) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon Q(a,y) = 0\}$
        Notation and definition from:
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.
        INPUT:
            'a' - Point in $\mathbb{P}^2$
        OUTPUT:
             A polynomial representing the fiber

        EXAMPLES::
           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
           sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
           sage: Y = x0*y0+x1*y1+x2*y2
           sage: X = WehlerK3Surface([Z,Y])
           sage: T = PP(1,1,0,1,0,0);
           sage: X.Qxa(T[0])
             5*y0^2 + 7*y0*y1 + y1^2 + 11*y1*y2 + y2^2
        """
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[1];
        PSYC = PSY.coordinate_ring()
        p = ASC.hom( [a[0],a[1],a[2]] + list( PSY.gens()), PSYC)
        return((p(self.Q)))
    
    def Sxa(self,a):
        r"""
        Function will return fiber defined by:
        $S^{x}_{a} = L^{x}_{a} \cap Q^{x}_{a}$
        Notation and definition from:
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.

        INPUT:
            'a' - Point in $\mathbb{P}^2$
              OUTPUT: A subscheme representing the fiber

        EXAMPLES::
           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
           sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
           sage: Y = x0*y0+x1*y1+x2*y2
           sage: X = WehlerK3Surface([Z,Y])
           sage: T = PP(1,1,0,1,0,0);
           sage: X.Sxa(T[0])
             Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              y0 + y1,
              5*y0^2 + 7*y0*y1 + y1^2 + 11*y1*y2 + y2^2
        """
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[1];
        PSYC = PSY.coordinate_ring()
        return PSY.subscheme([self.Lxa(a),self.Qxa(a)])

    def Lyb(self,b):
        r"""
        Function will return a fiber defined by:
        $L^{y}_{b} = \{(x,b) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon L(x,b) = 0\}$
        Notation and definition from:
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.
        INPUT:
            'a'- Point in Projective Space
        OUTPUT: 
            A polynomial representing the fiber
        EXAMPLES::
           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
           sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
           sage: Y = x0*y0+x1*y1+x2*y2
           sage: X = WehlerK3Surface([Z,Y])
           sage: T = PP(1,1,0,1,0,0);
           sage: X.Lyb(T[0])
             x0 + x1
        """
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[0];
        PSYC = PSY.coordinate_ring()
        p = ASC.hom(list(PSY.gens()) + [b[0],b[1],b[2]], PSYC)
        return ((p(self.L)))

    def Qyb(self,b):
        r"""
        Function will return a fiber defined by:
        $Q^{y}_{b} = \{(x,b) \in \mathbb{P}^{2} \times \mathbb{P}^{2} \colon Q(x,b) = 0\}$
        Notation and definition from:
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.
        INPUT:
            'b' - Point in Projective Space
        OUTPUT:
         A polynomial representing the fiber
        EXAMPLES::
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: T = PP(1,1,0,1,0,0);
            sage: X.Qyb(T[0])
             4*x0^2 + 6*x0*x1 + 3*x1^2 - x0*x2 - 4*x1*x2 - 2*x2^2
        """
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[0];
        PSYC = PSY.coordinate_ring()
        p = ASC.hom(list(PSY.gens()) + [b[0],b[1],b[2]], PSYC)
        return (p(self.Q))

    def Syb(self,b):
        r""" Function will return fiber defined by:
        $S^{y}_{b} = L^{y}_{b} \cap Q^{y}_{b}$
              Notation and definition from:
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.

        INPUT:
            'a' - Point in $\mathbb{P}^2$
              OUTPUT: A subscheme representing the fiber

        EXAMPLES::
           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
           sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
           sage: Y = x0*y0+x1*y1+x2*y2
           sage: X = WehlerK3Surface([Z,Y])
           sage: T = PP(1,1,0,1,0,0);
           sage: X.Syb(T[0])
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x0 + x1,
              4*x0^2 + 6*x0*x1 + 3*x1^2 - x0*x2 - 4*x1*x2 - 2*x2^2
        """
        AS = self.ambient_space()
        ASC = AS.coordinate_ring()
        PSY = AS[0];
        PSYC = PSY.coordinate_ring()
        return PSY.subscheme([self.Lyb(b),self.Qyb(b)])

    def Ramification_poly(self,i):
        r""" Function will return the Ramification polynomial defined by:
        $g^* = 2L^*_iL^*_jQ^*_{kk}$

        INPUT:
            'i' - Integer, either 0 (polynomial in y) or 1 (polynomial in x)

        OUTPUT:
            The Ramification polynomial for the Wehler K3 Surface
        
        EXAMPLES::
           sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
           sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
           sage: Y = x0*y0+x1*y1+x2*y2
           sage: X = WehlerK3Surface([Z,Y])
           sage: X.Ramification_poly(0)
           8*y0^5*y1 - 24*y0^4*y1^2 + 48*y0^2*y1^4 - 16*y0*y1^5 + y1^6 +\
            84*y0^3*y1^2*y2 + 46*y0^2*y1^3*y2 - 20*y0*y1^4*y2 + 16*y1^5*y2 +\
            53*y0^4*y2^2 + 56*y0^3*y1*y2^2 - 32*y0^2*y1^2*y2^2 - 80*y0*y1^3*y2^2 -\
            92*y1^4*y2^2 - 12*y0^2*y1*y2^3 - 168*y0*y1^2*y2^3 - 122*y1^3*y2^3 +\
            14*y0^2*y2^4 + 8*y0*y1*y2^4 - 112*y1^2*y2^4 + y2^6

        """
        return ((self._Lcoeff(i,0))**2)*(self._Qcoeff(i,1,2))**2+((self._Lcoeff(i,1))**2)*(self._Qcoeff(i,0,2)**2)+((self._Lcoeff(i,2))**2)*(self._Qcoeff(i,0,1)**2)-2*(self._Lcoeff(i,0))*(self._Lcoeff(i,1))*(self._Qcoeff(i,0,2))*(self._Qcoeff(i,1,2))-2*(self._Lcoeff(i,0))*(self._Lcoeff(i,2))*(self._Qcoeff(i,0,1))*(self._Qcoeff(i,1,2))-2*(self._Lcoeff(i,1))*(self._Lcoeff(i,2))*(self._Qcoeff(i,0,1))*(self._Qcoeff(i,0,2))+4*(self._Lcoeff(i,0))*(self._Lcoeff(i,1))*(self._Qcoeff(i,0,1))*(self._Qcoeff(i,2,2))+4*(self._Lcoeff(i,0))*(self._Lcoeff(i,2))*(self._Qcoeff(i,0,2))*(self._Qcoeff(i,1,1))+4*(self._Lcoeff(i,1))*(self._Lcoeff(i,2))*(self._Qcoeff(i,1,2))*(self._Qcoeff(i,0,0))-4*((self._Lcoeff(i,0))**2)*(self._Qcoeff(i,1,1))*(self._Qcoeff(i,2,2))-4*((self._Lcoeff(i,1))**2)*(self._Qcoeff(i,0,0))*(self._Qcoeff(i,2,2))-4*((self._Lcoeff(i,2))**2)*(self._Qcoeff(i,1,1))*(self._Qcoeff(i,0,0))

    @cached_method
    def is_degenerate(self):
        """
        Function will return True if there is a fiber (over the algebraic closure of the
        base ring) of dimension greater than 0 and False otherwise.
              INPUT: None.

        OUTPUT: Boolean value of True or False
        
        EXAMPLES::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.is_degenerate()
            True

          ::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.is_degenerate()
            False
              ::
                  #TODO: Should this really be true when the fiber is in an extension?
            #The problem is there is a degenerate fiber in an extension
            #PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],GF(3^10,'t'))
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],GF(3))
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.is_degenerate()
        """
        PP = self.ambient_space()
        K = FractionField(PP[0].base_ring())
        R=PP.coordinate_ring()
        PS = PP[0]; #check for x fibers
        vars = list(PS.gens())
        R0 = PolynomialRing(K,3,vars)
        I = R.ideal(self.Gpoly(1,0),self.Gpoly(1,1),self.Gpoly(1,2),self.Hpoly(1,0,1),self.Hpoly(1,0,2),self.Hpoly(1,1,2))
        phi = R.hom(vars+[0,0,0],R0)
        I = phi(I)
        if I.dimension() !=  0:
            return True

        PS = PP[1]; #check for y fibers
        vars = list(PS.gens())
        R0 = PolynomialRing(K,3,vars)
        I = R.ideal(self.Gpoly(0,0),self.Gpoly(0,1),self.Gpoly(0,2),self.Hpoly(0,0,1),self.Hpoly(0,0,2),self.Hpoly(0,1,2))
        phi = R.hom([0,0,0]+vars,R0)
        I = phi(I)
        if I.dimension() !=  0:
            return True
        return False

    def degenerate_fibers(self):
        r"""
        Function will return the (rational) degenerate fibers of the surface defined over the base ring,
        or the fraction field of the base ring if it is not a field.
              ALGORITHM:
              The criteria for degeneracy by the common vanishing of the polynomials $self.Gpoly(1,0)$, $self.Gpoly(1,1)$, $self.Gpoly(1,2)$, $self.Hpoly(1,0,1)$, $self.Hpoly(1,0,2)$, $self.Hpoly(1,1,2)$
        (for the first component) is from Proposition 1.4 in the following article:
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.
              This function finds the common solution through elimination via Groebner bases by using the .variety()
        function on the three affine charts in each component.

        INPUT: None

        OUTPUT:

        - The output is a list of lists where in the elements of lists are points in the appropriate projective space.
        The first list is the points whose pullback by the projection to the first component (projective space) is
        dimension greater than 0. The second list is points in the second component.
              NOTES::
              - Trac 13903 is needed for functionality over p-adic fields
              EXAMPLES::
                      sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
                sage: Y = x0*y0+x1*y1-x2*y2
                sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
                sage: X = WehlerK3Surface([Z,Y])
                sage: X.degenerate_fibers()
                [[], [(1 : 0 : 0)]]
              ::
                      sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
                sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
                sage: Y = x0*y0+x1*y1+x2*y2
                sage: X = WehlerK3Surface([Z,Y])
                sage: X.degenerate_fibers()
                [[], []]
                      ::
                      sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
                sage: R = PP.coordinate_ring()
                sage: l = y0*x0 + y1*x1 + (y0 - y1)*x2
                sage: q = (y1*y0 + y2^2)*x0^2 + ((y0^2 - y2*y1)*x1 + (y0^2 + (y1^2 - y2^2))*x2)*x0 + (y2*y0 + y1^2)*x1^2+ (y0^2 + (-y1^2 + y2^2))*x2*x1
                sage: X = WehlerK3Surface([l,q])
                sage: X.degenerate_fibers()
                [[(-1 : 1 : 1), (0 : 0 : 1)], [(-1 : -1 : 1), (0 : 0 : 1)]]
        """
        PP = self.ambient_space()
        PSX = PP[0];
        vars = list(PSX.gens())
        K = FractionField(PSX.base_ring())
        R0 = PolynomialRing(K,3,vars)
        I = ideal(self.Gpoly(1,0),self.Gpoly(1,1),self.Gpoly(1,2),self.Hpoly(1,0,1),self.Hpoly(1,0,2),self.Hpoly(1,1,2))
        phi = PP.coordinate_ring().hom(vars+[0,0,0],R0)
        I = phi(I)
        xFibers = []
        #check affine charts
        for n in range(3):
            affvars = list(R0.gens())
            del affvars[n]
            R1 = PolynomialRing(K,2,affvars,order = 'lex')
            mapvars = list(R1.gens())
            mapvars.insert(n,1)
            phi1 = R0.hom(mapvars,R1)
            J = phi1(I)
            if (J.dimension() == 0):
                Var = J.variety()
                if Var !=  [{}]:
                    for d in Var: #iterate through dictionaries
                        P = [] #new point
                        for z in mapvars: #assign coordinate values
                            if (z == 1):
                                P.append(1)
                            else:
                                P.append(d[z])
                        MP = PSX(P) #make projective point
                        if (MP in xFibers) == False:
                            xFibers.append(MP)
        PSY = PP[1]
        vars = list(PSY.gens())
        K = FractionField(PSY.base_ring())
        R0 = PolynomialRing(K,3,vars)
        I = ideal(self.Gpoly(0,0),self.Gpoly(0,1),self.Gpoly(0,2),self.Hpoly(0,0,1),self.Hpoly(0,0,2),self.Hpoly(0,1,2))
        phi = PP.coordinate_ring().hom([0,0,0]+vars,R0)
        I = phi(I)
        yFibers = []
        #check affine charts
        for n in range(3):
            affvars = list(R0.gens())
            del affvars[n]
            R1 = PolynomialRing(K,2,affvars,order = 'lex')
            mapvars = list(R1.gens())
            mapvars.insert(n,1)
            phi1 = R0.hom(mapvars,R1)
            J = phi1(I)
            if (J.dimension() == 0):
                Var = J.variety()
                if Var !=  [{}]:
                    for d in Var: #iterate through dictionaries
                        P = [] #new point                         for z in mapvars: #assign coordinate values
                        if (z == 1):
                            P.append(1)
                        else:
                            P.append(d[z])
                        MP = PSY(P) #make projective point
                        if (MP in yFibers) == False:
                            yFibers.append(MP)
        return [xFibers,yFibers]

    @cached_method
    def degenerate_primes(self,check = True):
        r"""
        Determine which primes $p$ self has degenerate fibers over $GF(p)$. If check is False, then
        may return primes that do not have degenerate fibers. Raises an error if the surface is degenerate.
        
        INPUT:
        
        - check Boolean - True then the primes are verified
        
        OUTPUT:
        
        - list of primes.
        
        EXAMPLES::
        
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ,6)
            sage: L =  y0*x0 + (y1*x1 + y2*x2)
            sage: Q = (2*y0^2 + y2*y0 + (2*y1^2 + y2^2))*x0^2 + ((y0^2 + y1*y0 + (y1^2 + 2*y2*y1 + y2^2))*x1 + (2*y1^2 + y2*y1 + y2^2)*x2)*x0 + ((2*y0^2 + (y1 + 2*y2)*y0 + (2*y1^2 + y2*y1))*x1^2 + ((2*y1 + 2*y2)*y0 + (y1^2 + y2*y1 + 2*y2^2))*x2*x1 + (2*y0^2 + y1*y0 + (2*y1^2 + y2^2))*x2^2)
            sage: X = WehlerK3Surface([L,Q])
            sage: X.degenerate_primes()
            [2, 3, 5, 11, 23, 47, 48747691, 111301831]
        """

        if self.ambient_space().base_ring() !=  ZZ and self.ambient_space().base_ring() !=  QQ:
            raise TypeError, "Must be ZZ or QQ"
        if self.is_degenerate():
            raise TypeError, "Surface is degenerate at all primes"
        PP = self.ambient_space()

        #x-fibers
        PSX = PP[0];
        vars = list(PSX.gens())
        K = PSX.base_ring()
        R = PolynomialRing(K,3,vars)
        I = ideal(self.Gpoly(1,0),self.Gpoly(1,1),self.Gpoly(1,2),self.Hpoly(1,0,1),self.Hpoly(1,0,2),self.Hpoly(1,1,2))
        phi = PP.coordinate_ring().hom(vars+[0,0,0],R)
        I = phi(I)
        #if J.dimension()>0:
        #    raise TypeError, "Not a morphism."
        badPrimes = []

        #move the ideal to the ring of integers
        if R.base_ring().is_field():
            S = PolynomialRing(R.base_ring().ring_of_integers(),R.gens(),R.ngens())
            I = S.ideal(I.gens())
        GB = I.groebner_basis()
        #get the primes dividing the coefficients of the monomials x_i^k_i
        for i in range(len(GB)):
            LT = GB[i].lt().degrees()
            power = 0
            for j in range(R.ngens()):
                if LT[j] != 0:
                    power += 1
            if power == 1:
                badPrimes = badPrimes+GB[i].lt().coefficients()[0].support()

        #y-fibers
        PSY = PP[1];
        vars = list(PSY.gens())
        K = PSY.base_ring()
        R = PolynomialRing(K,3,vars)
        I = ideal(self.Gpoly(0,0),self.Gpoly(0,1),self.Gpoly(0,2),self.Hpoly(0,0,1),self.Hpoly(0,0,2),self.Hpoly(0,1,2))
        phi = PP.coordinate_ring().hom([0,0,0]+vars,R)
        I = phi(I)
        #move the ideal to the ring of integers
        if R.base_ring().is_field():
            S = PolynomialRing(R.base_ring().ring_of_integers(),R.gens(),R.ngens())
            I = S.ideal(I.gens())
        GB = I.groebner_basis()
        #get the primes dividing the coefficients of the monomials x_i^k_i
        for i in range(len(GB)):
            LT = GB[i].lt().degrees()
            power = 0
            for j in range(R.ngens()):
                if LT[j] != 0:
                    power += 1
            if power == 1:
                badPrimes = badPrimes+GB[i].lt().coefficients()[0].support()
        badPrimes = list(set(badPrimes))
        badPrimes.sort()
        #check to return only the truly bad primes
        if check == True:
            for p in badPrimes:
                X = self.change_ring(GF(p))
                if X.is_degenerate() == False:
                    badPrimes.remove(p)
        return(badPrimes)
    def is_smooth(self):
        r"""
        Function will return the status of the smoothness of the surface
              ALGORITHM:
              The checks to confirm that none of the 2x2 minors of the jacobian generated from the Biquadratic and Bilinear forms have
        no common vanishing points
              INPUT: None.
              OUTPUT:
              - Boolean
              EXAMPLES::
                  sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = x0^2*y0*y1 + x0^2*y2^2 - x0*x1*y1*y2 + x1^2*y2*y1 +x2^2*y2^2 + x2^2*y1^2 +x1^2*y2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.is_smooth()
            False
                  ::
                      sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            True"""

        vars = list(self.ambient_space().gens())
        M = jacobian([self.L,self.Q], vars)
        R = self.ambient_space().coordinate_ring()
        I = R.ideal(M.minors(2)+[self.L,self.Q])
        T = PolynomialRing(self.ambient_space().base_ring().fraction_field(),4, 'h')
        #check the 9 affine charts for a singular point
        for l in mrange([3,3]):
            vars = list(T.gens())
            vars.insert(l[0],1)
            vars.insert(3+l[1],1)
            phi = R.hom(vars,T)
            J = phi(I)
            if J.dimension() != -1:
                return False
        return True

    def sigmaX(self,P, check = True):
        r""" Function returns the involution on the surface S induced by the double covers

        ALGORITHM:

        Refer to Section 6: "An algorithm to compute $\sigma_x$, $\sigma_y$, $\phi$, and $\psi$" in
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.

        INPUT:
            'P' - a point in $\mathbb{P}^2 \times $\mathbb{P}^2$
        OUTPUt:
            A point in $\mathbb{P}^2 \times $\mathbb{P}^2$
              EXAMPLES::
                  sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: T = PP(0,0,1,1,0,0)
            sage: X.sigmaX(T)
            (0 : 0 : 1 , 0 : 2 : 0)

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: l = y0*x0 + y1*x1 + (y0 - y1)*x2
            sage: q = (y1*y0)*x0^2 + ((y0^2)*x1 + (y0^2 + (y1^2 - y2^2))*x2)*x0 + (y2*y0 + y1^2)*x1^2+ (y0^2 + (-y1^2 + y2^2))*x2*x1
            sage: X = WehlerK3Surface([l,q])
            sage: X.sigmaX(X([1,0,0,0,1,-2]))
            (1 : 0 : 0 , 0 : -16 : -32)
            sage: X.sigmaX(X([1,0,0,0,0,1])) #s = infty
            (1 : 0 : 0 , 0 : 0 : -1)
            sage: X.sigmaX(X([-1,1,1,-1,-1,1])) #s = 0
            (-1 : 1 : 1 , -2 : -2 : -1)
            sage: X.sigmaX(X([0,0,1,1,1,0])) #s = infty
            (0 : 0 : 1 , 1 : 1 : 0)
            sage: X.sigmaX(X([0,0,1,1,1,1]))
            (0 : 0 : 1 , 2 : 2 : -2)
        """
        if check:
            self._check_satisfies_equations(P)
        pt = list(P[0]) + [0,0,0]
        if(P[1][0] != 0):
            [a,b,c] = [P[1][0]*self.Gpoly(1,0)(*pt),\
                       -1*P[1][0]*self.Hpoly(1,0,1)(*pt)-P[1][1]*self.Gpoly(1,0)(*pt),\
                       -P[1][0]*self.Hpoly(1,0,2)(*pt)-P[1][2]*self.Gpoly(1,0)(*pt)]
        elif(P[1][1] != 0):
            [a,b,c] = [-1*P[1][1]*self.Hpoly(1,0,1)(*pt)-P[1][0]*self.Gpoly(1,1)(*pt),\
                        P[1][1]*self.Gpoly(1,1)(*pt),\
                       -P[1][1]*self.Hpoly(1,1,2)(*pt)-P[1][2]*self.Gpoly(1,1)(*pt)]
        elif(P[1][2]!= 0):
            [a,b,c] = [-1*P[1][2]*self.Hpoly(1,0,2)(*pt)-P[1][0]*self.Gpoly(1,2)(*pt),\
                       -P[1][2]*self.Hpoly(1,1,2)(*pt)-P[1][1]*self.Gpoly(1,2)(*pt),\
                       P[1][2]*self.Gpoly(1,2)(*pt)]
        Point = [P[0][0],P[0][1],P[0][2],a,b,c]
        if [a,b,c]!= [0,0,0]:
            return self.point(Point,False)
        R = self.ambient_space().coordinate_ring()
        BR = self.ambient_space().base_ring()
        if P[0][0]!= 0:
            S = PolynomialRing(BR,6,'z')
            s0,s1,w1,z0,z1,z2 = S.gens()
            t = w1-BR(P[0][1]/P[0][0])
            t1 = BR(P[0][1]/P[0][0])
            phi = R.hom([s0,s0*w1,s1*t + s0*P[0][2]/P[0][0],z0,z1,z2],S)
            T = [phi(self.L),phi(self.Q),\
                 phi(self.Gpoly(1,0)),\
                 phi(self.Gpoly(1,1)),\
                 phi(self.Gpoly(1,2)),\
                -phi(self.Hpoly(1,0,1)),\
                -phi(self.Hpoly(1,0,2)),\
                -phi(self.Hpoly(1,1,2))]
        elif P[0][1]!= 0:
            S = PolynomialRing(BR,6,'z')
            s0,s1,w1,z0,z1,z2 = S.gens()
            t = w1-BR(P[0][0]/P[0][1])
            t1 = BR(P[0][0]/P[0][1])
            phi = R.hom([s0*w1,s0,s1*t + s0*P[0][2]/P[0][1],z0,z1,z2],S)
            T = [phi(self.L),\
                 phi(self.Q),\
                 phi(self.Gpoly(1,0)),\
                 phi(self.Gpoly(1,1)),\
                 phi(self.Gpoly(1,2)),\
                -phi(self.Hpoly(1,0,1)),\
                -phi(self.Hpoly(1,0,2)),\
                -phi(self.Hpoly(1,1,2))]
        elif P[0][2]!= 0:
            S = PolynomialRing(BR,6,'z')
            s0,s1,w1,z0,z1,z2 = S.gens()
            t = w1-BR(P[0][1]/P[0][2])
            t1 = BR(P[0][1]/P[0][2])
            phi = R.hom([s1*(t) + s0*P[0][0]/P[0][2],s0*w1,s0,z0,z1,z2],S)
            T = [phi(self.L),\
                 phi(self.Q),\
                 phi(self.Gpoly(1,0)),\
                 phi(self.Gpoly(1,1)),\
                 phi(self.Gpoly(1,2)),\
                -phi(self.Hpoly(1,0,1)),\
                -phi(self.Hpoly(1,0,2)),\
                -phi(self.Hpoly(1,1,2))]
        maxexp = []

        for i in range(2,len(T)):
            e = 0
            while (T[i]/t^e).subs({w1:t1}) == 0:
                e+= 1
            maxexp.append(e)
            e  = min(maxexp)
        for i in range(0,2):
            while T[i].subs({w1:t1}) == 0:
                T[i] = T[i]/t
            T[i] = T[i].subs({w1:t1})
        for i in range(2,len(T)):
            T[i] = T[i]/t^e
            T[i] = T[i].subs({w1:t1})

        RR = PolynomialRing(BR,5,'z',order = 'lex')
        s0,s1,z0,z1,z2 = RR.gens()
        I = RR.ideal([RR(T[0]),\
                      RR(T[1]),\
                      RR(T[2]) - P[1][0]*z0, RR(T[3])-P[1][1]*z1,RR(T[4])-P[1][2]*z2,\
                      RR(T[5]) - (P[1][0]*z1 + P[1][1]*z0),\
                      RR(T[6]) - (P[1][0]*z2 + P[1][2]*z0),\
                      RR(T[7]) - (P[1][1]*z2 + P[1][2]*z1)])

        SS = PolynomialRing(BR,4,'z',order = 'lex')
        s,z0,z1,z2 = SS.gens()
        phi = RR.hom([s,1,z0,z1,z2],SS)
        V = phi(I).variety()
        if len(V) !=  1:
            for D in V:
                if D[s]!= 0: #(0,0,0,0) is always a solution
                    [a,b,c] = [D[z0],D[z1],D[z2]]
        else: #maybe s == 0, but we need to kill the identically 0 portions
            newT = [phi(t) for t in T]
            for i in range(2):
                while newT[i]!= 0 and s.divides(newT[i]):
                    newT[i] = SS(newT[i]/s)
            maxexp = []
            for i in range(2,len(T)):
                e = 0
                if newT[i]!= 0:
                    while (newT[i]/s^e).subs({s:0}) == 0:
                        e+= 1
                    maxexp.append(e)
            e = min(maxexp)
            for i in range(2,len(T)):
                newT[i] = newT[i]/s^e
            II = SS.ideal([SS(newT[0]),\
                           SS(newT[1]),\
                           SS(newT[2]) - P[1][0]*z0,\
                           SS(newT[3])-P[1][1]*z1,\
                           SS(newT[4])-P[1][2]*z2,\
                           SS(newT[5]) - (P[1][0]*z1 + P[1][1]*z0),\
                           SS(newT[6]) - (P[1][0]*z2 + P[1][2]*z0),\
                           SS(newT[7]) - (P[1][1]*z2 + P[1][2]*z1)])
            SSS = PolynomialRing(BR,3,'z',order = 'lex')
            z0,z1,z2 = SSS.gens()
            phi = SS.hom([0,z0,z1,z2],SSS)
            V = phi(II).variety()
            if len(V) != 0:
                [a,b,c] = [V[0][z0],V[0][z1],V[0][z2]]
            if len(V) == 0 or [a,b,c] == [0,0,0]:
                SS = PolynomialRing(BR,3,'z',order = 'lex')
                z0,z1,z2 = SS.gens()
                phi = RR.hom([1,0,z0,z1,z2],SS)
                V = phi(I).variety()
                [a,b,c] = [V[0][z0],V[0][z1],V[0][z2]]
        Point = [P[0][0],P[0][1],P[0][2],a,b,c]
        return self.point(Point,False)

    def sigmaY(self,P, check = True):
        r""" Function returns the involution on the surface S induced by the double covers

        ALGORITHM:

        Refer to Section 6: "An algorithm to compute $\sigma_x$, $\sigma_y$, $\phi$, and $\psi$" in
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.

        INPUT:
            'P' - a point in $\mathbb{P}^2 \times $\mathbb{P}^2$
        OUTPUt:
            A point in $\mathbb{P}^2 \times $\mathbb{P}^2$
              EXAMPLES::
                      sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: T = PP(0,0,1,1,0,0)
            sage: X.sigmaY(T)
            (0 : 0 : 1 : 1 : 0 : 0)
                  degenerate examples
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: l = y0*x0 + y1*x1 + (y0 - y1)*x2
            sage: q = (y1*y0)*x0^2 + ((y0^2)*x1 + (y0^2 + (y1^2 - y2^2))*x2)*x0 + (y2*y0 + y1^2)*x1^2+ (y0^2 + (-y1^2 + y2^2))*x2*x1
            sage: X = WehlerK3Surface([l,q])
            sage: X.sigmaY(X([1,-1,0,-1,-1,1]))
            (-1/256 : 1/256 : -5/128 : -1 : -1 : 1)
            sage: X.sigmaY(X([0,0,1,-1,-1,1]))
            (-4 : 4 : 1 : -1 : -1 : 1)
            sage: X.sigmaY(X([1,2,0,0,0,1]))
            (3/8 : 3/8 : -1/8 : 0 : 0 : 1)
            sage: X.sigmaY(X([1,1,1,0,0,1])) #s = infty
            (1 : 0 : 0 : 0 : 0 : 1)
        """
        if check:
            self._check_satisfies_equations(P)
        pt = [0,0,0] + list(P[1])
        if(P[0][0] != 0):
            [a,b,c] = [P[0][0]*self.Gpoly(0,0)(*pt),\
                      -1*P[0][0]*self.Hpoly(0,0,1)(*pt) - P[0][1]*self.Gpoly(0,0)(*pt),\
                      -P[0][0]*self.Hpoly(0,0,2)(*pt) - P[0][2]*self.Gpoly(0,0)(*pt)]
        elif(P[0][1] != 0):
            [a,b,c] = [-1*P[0][1]*self.Hpoly(0,0,1)(*pt) - P[0][0]*self.Gpoly(0,1)(*pt),\
                       P[0][1]*self.Gpoly(0,1)(*pt),\
                       -P[0][1]*self.Hpoly(0,1,2)(*pt) - P[0][2]*self.Gpoly(0,1)(*pt)]
        elif(P[0][2] != 0):
            [a,b,c] = [-1*P[0][2]*self.Hpoly(0,0,2)(*pt) - P[0][0]*self.Gpoly(0,2)(*pt),\
                       -P[0][2]*self.Hpoly(0,1,2)(*pt) - P[0][1]*self.Gpoly(0,2)(*pt),\
                       P[0][2]*self.Gpoly(0,2)(*pt)]
        Point = [a,b,c,P[1][0],P[1][1],P[1][2]]
        if [a,b,c]!= [0,0,0]:
            return self.point(Point,False)
        R = self.ambient_space().coordinate_ring()
        BR = self.ambient_space().base_ring()
        if P[1][0]!= 0:
            S = PolynomialRing(BR,6,'z')
            z0,z1,z2,s0,s1,w1 = S.gens()
            t = w1-BR(P[1][1]/P[1][0])
            t1 = BR(P[1][1]/P[1][0])
            phi = R.hom([z0,z1,z2,s0,s0*w1,s1*t + s0*P[1][2]/P[1][0]],S)
            T = [phi(self.L),\
                 phi(self.Q),\
                 phi(self.Gpoly(0,0)),\
                 phi(self.Gpoly(0,1)),\
                 phi(self.Gpoly(0,2)),\
                -phi(self.Hpoly(0,0,1)),\
                -phi(self.Hpoly(0,0,2)),\
                -phi(self.Hpoly(0,1,2))]
        elif P[1][1]!= 0:
            S = PolynomialRing(BR,6,'z')
            z0,z1,z2,s0,s1,w1 = S.gens()
            t = w1-BR(P[1][0]/P[1][1])
            t1 = BR(P[1][0]/P[1][1])
            phi = R.hom([z0,z1,z2,s0*w1,s0,s1*t + s0*P[1][2]/P[1][1]],S)
            T = [phi(self.L),\
                 phi(self.Q),\
                 phi(self.Gpoly(0,0)),\
                 phi(self.Gpoly(0,1)),\
                 phi(self.Gpoly(0,2)),\
                -phi(self.Hpoly(0,0,1)),\
                -phi(self.Hpoly(0,0,2)),\
                -phi(self.Hpoly(0,1,2))]
        elif P[1][2]!= 0:
            S = PolynomialRing(BR,6,'z')
            z0,z1,z2,s0,s1,w1 = S.gens()
            t = w1-BR(P[1][1]/P[1][2])
            t1 = BR(P[1][1]/P[1][2])
            phi = R.hom([z0,z1,z2,s1*(t) + s0*P[1][0]/P[1][2],s0*w1,s0],S)
            T = [phi(self.L),\
                 phi(self.Q),\
                 phi(self.Gpoly(0,0)),\
                 phi(self.Gpoly(0,1)),\
                 phi(self.Gpoly(0,2)),\
                -phi(self.Hpoly(0,0,1)),\
                -phi(self.Hpoly(0,0,2)),\
                -phi(self.Hpoly(0,1,2))]
            maxexp = []

        for i in range(2,len(T)):
            e = 0
            while (T[i]/t^e).subs({w1:t1}) == 0:
                e+= 1
            maxexp.append(e)
        e = min(maxexp)

        for i in range(0,2):
            while T[i].subs({w1:t1}) == 0:
                T[i] = T[i]/t
            T[i] = T[i].subs({w1:t1})
        for i in range(2,len(T)):
            T[i] = T[i]/t^e
            T[i] = T[i].subs({w1:t1})

        RR = PolynomialRing(BR,5,'z',order = 'lex')
        s0,s1,z0,z1,z2 = RR.gens()
        I = RR.ideal([RR(T[0]),\
                      RR(T[1]),\
                      RR(T[2]) - P[0][0]*z0,\
                      RR(T[3])-P[0][1]*z1,\
                      RR(T[4])-P[0][2]*z2,\
                      RR(T[5]) - (P[0][0]*z1 + P[0][1]*z0),\
                      RR(T[6]) - (P[0][0]*z2 + P[0][2]*z0),\
                      RR(T[7]) - (P[0][1]*z2 + P[0][2]*z1)])
        SS = PolynomialRing(BR,4,'z',order = 'lex')
        s,z0,z1,z2 = SS.gens()
        phi = RR.hom([s,1,z0,z1,z2],SS)
        V = phi(I).variety()

        if len(V) !=  1:
            for D in V:
                if D[s]!= 0:
                    [a,b,c] = [D[z0],D[z1],D[z2]]
        else:
            newT = [phi(t) for t in T]
            for i in range(2):
                while newT[i]!= 0 and s.divides(newT[i]):
                    newT[i] = SS(newT[i]/s)
            maxexp = []
            for i in range(2,len(T)):
                e = 0
                if newT[i]!= 0:
                    while (newT[i]/s^e).subs({s:0}) == 0:
                        e+= 1
                    maxexp.append(e)
            e = min(maxexp)
            for i in range(2,len(T)):
                newT[i] = newT[i]/s^e
            II = SS.ideal([SS(newT[0]),\
                           SS(newT[1]),\
                           SS(newT[2]) - P[0][0]*z0,\
                            SS(newT[3])-P[0][1]*z1,\
                            SS(newT[4])-P[0][2]*z2,\
                            SS(newT[5]) - (P[0][0]*z1 + P[0][1]*z0),\
                            SS(newT[6]) - (P[0][0]*z2 + P[0][2]*z0),\
                            SS(newT[7]) - (P[0][1]*z2 + P[0][2]*z1)])
            SSS = PolynomialRing(BR,3,'z',order = 'lex')
            z0,z1,z2 = SSS.gens()
            phi = SS.hom([0,z0,z1,z2],SSS)
            V = phi(II).variety()

            if len(V) != 0:
                [a,b,c] = [V[0][z0],V[0][z1],V[0][z2]]
            if len(V) == 0 or [a,b,c] == [0,0,0]:
                SS = PolynomialRing(BR,3,'z',order = 'lex')
                z0,z1,z2 = SS.gens()
                phi = RR.hom([1,0,z0,z1,z2],SS)
                V = phi(I).variety()
                [a,b,c] = [V[0][z0],V[0][z1],V[0][z2]]

        Point = [a,b,c,P[1][0],P[1][1],P[1][2]]
        return self.point(Point,False)

    def phi(self,a,check = True ):
        r"""
        Evaluates the function $\phi = \sigma_y \circ \sigma_x$

        ALGORITHM:
          Refer to Section 6: "An algorithm to compute $\sigma_x$, $\sigma_y$, $\phi$, and $\psi$" in
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.

        INPUT:

        - ``a`` - Point in $\mathbb{P} \times \mathbb{P}$

        OUTPUT:

        Point in $\mathbb{P} \times \mathbb{P}$

        EXAMPLES:::

        sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
        sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
        sage: Y = x0*y0+x1*y1+x2*y2
        sage: X = WehlerK3Surface([Z,Y])
        sage: T = PP(0,0,1,1,0,0)
        sage: X.phi(T)             
        (16 : 0 : -16 : 0 : 2 : 0)


        """

        A = self.sigmaX(a, check)
        return self.sigmaY(A, check)

    def psi(self,a, check=True):
        r"""
        Evaluates the function $\psi = \sigma_x \circ \sigma_y$

        ALGORITHM:
          Refer to Section 6: "An algorithm to compute $\sigma_x$, $\sigma_y$, $\phi$, and $\psi$" in
        G. Call and J. Silverman. Computing the canonical height on K3 surfaces. Math. Comp., 65:(259-290), 1996.

        INPUT:

        - ``a`` - Point in $\mathbb{P} \times \mathbb{P}$

        OUTPUT:

        Point in $\mathbb{P} \times \mathbb{P}$

        EXAMPLES:::

        sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
        sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
        sage: Y = x0*y0+x1*y1+x2*y2
        sage: X = WehlerK3Surface([Z,Y])
        sage: T = PP(0,0,1,1,0,0)
        sage: X.psi(T)            (0 : 0 : 1 : 0 : 2 : 0)
        """

        A = self.sigmaY(a, check)
        return self.sigmaX(A, check)

    def lambda_plus(self,P,v,N,m,n, prec = 100):
        r"""
        Evaluates the `+` local canonical height function of Call-Silverman at the place ``v``
        for ``P`` with ``N`` terms of the series.
              Use ``v = 0`` for the archimedean place. Must be over `\ZZ` or `\QQ`.
              ALGORITHM:
              See ``Computing the Canonical Height on K3 Surfaces``, G. Call and J. Silverman,
        Mathematics of Computation, 65(213):259-290, 1996.
              INPUT:
              - ``P`` - a projective point
              - ``N`` - positive integer. number of terms of the series to use
              - ``v`` - non-negative integer. a place, use v = 0 for the archimedean place
              - ``m,n`` - positive integers, We compute the local heigt for the divisor $E_{mn}^{+}$.
                    These must be indices of non-zero coordinates of the point ``P``.
              - ``prec`` - float point or p-adic precision, default: 100
              OUTPUT:
              - a real number

        EXAMPLES:::
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: P = X([0,0,1,1,0,0])
            sage: X.lambda_plus(P,0,10,2,0)
            0.892307051691616
        """
        if not (v == 0 or is_prime(v)):
            raise ValueError("Invalid valuation ( = %s) entered."%v)
        R = RealField(prec)
        if v == 0:
            K = R
        else:
            K = Qp(v, prec)
        PK = P.change_ring(K)
        W = self.change_ring(K)
        Rx = W.ambient_space().coordinate_ring().hom(\
                                                     list(W.ambient_space()[0].coordinate_ring().gens())+[0,0,0],\
                                                     W.ambient_space()[0].coordinate_ring())
        Ry = W.ambient_space().coordinate_ring().hom(\
                                                     [0,0,0]+list(W.ambient_space()[1].coordinate_ring().gens()),\
                                                     W.ambient_space()[1].coordinate_ring())
        beta = R(2+sqrt(3))
        L = [x.abs() for x in list(PK[0])]
        i = L.index(max(L))
        L = [y.abs() for y in list(PK[1])]
        j = L.index(max(L))
        localHeight = beta*R((PK[0][i]/PK[0][m]).abs()).log() - R((PK[1][j]/PK[1][n]).abs()).log()
        for e in range(N):
            Q = W.phi(PK, check=False)
            L = [x.abs() for x in list(Q[0])]
            k = L.index(max(L))
            L = [y.abs() for y in list(Q[1])]
            l = L.index(max(L))
            newP = copy(PK)
            newP.scale_by([1/PK[0][i],1])
            if PK[1][j].abs() <= PK[1][l].abs():
                if l == 0:
                    B = Rx(W.Gpoly(1,0))(tuple(newP[0]))*PK[1][j]/PK[1][l]
                elif l == 1:
                    B = Rx(W.Gpoly(1,1))(tuple(newP[0]))*PK[1][j]/PK[1][l]
                else:
                    B = Rx(W.Gpoly(1,2))(tuple(newP[0]))*PK[1][j]/PK[1][l]
            else:
                if j == 0:
                    B = -Rx(W.Gpoly(1,0))(tuple(newP[0]))*PK[1][l]/PK[1][j]
                    if l == 1:
                        B = B-Rx(W.Hpoly(1,0,1))(tuple(newP[0]))
                    else:
                        B = B-Rx(W.Hpoly(1,0,2))(tuple(newP[0]))
                elif j == 1:
                    B = -Rx(W.Gpoly(1,1))(tuple(newP[0]))*PK[1][l]/PK[1][j]
                    if l == 0:
                        B = B-Rx(W.Hpoly(1,0,1))(tuple(newP[0]))
                    else:
                        B = B-Rx(W.Hpoly(1,1,2))(tuple(newP[0]))
                else:
                    B = -Rx(W.Gpoly(1,2))(tuple(newP[0]))*PK[1][l]/PK[1][j]
                    if l == 0:
                        B = B-Rx(W.Hpoly(1,0,2))(tuple(newP[0]))
                    else:
                        B = B-Rx(W.Hpoly(1,1,2))(tuple(newP[0]))
            newQ = copy(Q)
            newQ.scale_by([1,1/Q[1][l]])
            if PK[0][i].abs() <= PK[0][k].abs():
                if k == 0:
                    A = Ry(W.Gpoly(0,0))(tuple(newQ[1]))*Q[0][i]/Q[0][k]
                elif k == 1:
                    A = Ry(W.Gpoly(0,1))(tuple(newQ[1]))*Q[0][i]/Q[0][k]
                else:
                    A = Ry(W.Gpoly(0,2))(tuple(newQ[1]))*Q[0][i]/Q[0][k]
            else:
                if i == 0:
                    A = -Ry(W.Gpoly(0,0))(tuple(newQ[1]))*PK[0][k]/PK[0][i]
                    if k == 1:
                        A = A-Ry(W.Hpoly(0,0,1))(tuple(newQ[1]))
                    else:
                        A = A-Ry(W.Hpoly(0,0,2))(tuple(newQ[1]))
                elif i == 1:
                    A = -Ry(W.Gpoly(0,1))(tuple(newQ[1]))*PK[0][k]/PK[0][i]
                    if k == 0:
                        A = A-Ry(W.Hpoly(0,0,1))(tuple(newQ[1]))
                    else:
                        A = A-Ry(W.Hpoly(0,1,2))(tuple(newQ[1]))
                else:
                    A = -Ry(W.Gpoly(0,2))(tuple(newQ[1]))*PK[0][k]/PK[0][i]
                    if k == 0:
                        A = A-Ry(W.Hpoly(0,0,2))(tuple(newQ[1]))
                    else:
                        A = A-Ry(W.Hpoly(0,1,2))(tuple(newQ[1]))
            localHeight += beta**(-2*R(e)-1)*R(A.abs()).log() + beta**(-2*R(e))*R(B.abs()).log()
            i = k
            j = l
            newQ.scale_by([1/Q[0][k],1])
            PK = newQ
        return localHeight

    def lambda_minus(self,P,v,N,m,n,prec = 100):
        r"""
        Evaluates the `-` local canonical height function of Call-Silverman at the place ``v``
        for ``P`` with ``N`` terms of the series.
              Use ``v = 0`` for the archimedean place. Must be over `\ZZ` or `\QQ`.
              ALGORITHM:
              See ``Computing the Canonical Height on K3 Surfaces``, G. Call and J. Silverman,
        Mathematics of Computation, 65(213):259-290, 1996.
              INPUT:
              - ``P`` - a projective point
              - ``N`` - positive integer. number of terms of the series to use
              - ``v`` - non-negative integer. a place, use v = 0 for the archimedean place
              - ``m,n`` - positive integers, We compute the local heigt for the divisor $E_{mn}^{+}$.
                    These must be indices of non-zero coordinates of the point ``P``.
              - ``prec`` - float point or p-adic precision, default: 100
              OUTPUT:
              - a real number

        EXAMPLES:::
                      sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: P = X([0,0,1,1,0,0])
            sage: X.lambda_minus(P,2,20,2,0,200)
            -0.18573351672047135037172805779671791488351056677474271893705
        """
        R = RealField(prec)
        if v == 0:
            K = R
        else:
            K = Qp(v, prec)
        PK = P.change_ring(K)
        W = self.change_ring(K)
        Rx = W.ambient_space().coordinate_ring().hom(list(W.ambient_space()[0].coordinate_ring().gens())+[0,0,0],W.ambient_space()[0].coordinate_ring())
        Ry = W.ambient_space().coordinate_ring().hom([0,0,0]+list(W.ambient_space()[1].coordinate_ring().gens()),W.ambient_space()[1].coordinate_ring())
        beta = R(2+sqrt(3))
        L = [x.abs() for x in list(PK[0])]
        j = L.index(max(L))
        L = [y.abs() for y in list(PK[1])]
        i = L.index(max(L))
        localHeight = beta*R((PK[1][i]/PK[1][n]).abs()).log() - R((PK[0][j]/PK[0][m]).abs()).log()
        for e in range(N):
            Q = W.psi(PK, check=False)
            L = [x.abs() for x in list(Q[0])]
            l = L.index(max(L))
            L = [y.abs() for y in list(Q[1])]
            k = L.index(max(L))
            newP = copy(PK)
            newP.scale_by([1,1/PK[1][i]])
            if PK[0][j].abs() <= PK[0][l].abs():
                if l == 0:
                    B = Ry(W.Gpoly(0,0))(tuple(newP[1]))*PK[0][j]/PK[0][l]
                elif l == 1:
                    B = Ry(W.Gpoly(0,1))(tuple(newP[1]))*PK[0][j]/PK[0][l]
                else:
                    B = Ry(W.Gpoly(0,2))(tuple(newP[1]))*PK[0][j]/PK[0][l]
            else:
                if j == 0:
                    B = -Ry(W.Gpoly(0,0))(tuple(newP[1]))*PK[0][l]/PK[0][j]
                    if l == 1:
                        B = B-Ry(W.Hpoly(0,0,1))(tuple(newP[1]))
                    else:
                        B = B-Ry(W.Hpoly(0,0,2))(tuple(newP[1]))
                elif j == 1:
                    B = -Ry(W.Gpoly(0,1))(tuple(newP[1]))*PK[0][l]/PK[0][j]
                    if l == 0:
                        B = B-Ry(W.Hpoly(0,0,1))(tuple(newP[1]))
                    else:
                        B = B-Ry(W.Hpoly(0,1,2))(tuple(newP[1]))
                else:
                    B = -Ry(W.Gpoly(0,2))(tuple(newP[1]))*PK[0][l]/PK[0][j]
                    if l == 0:
                        B = B-Ry(W.Hpoly(0,0,2))(tuple(newP[1]))
                    else:
                        B = B-Ry(W.Hpoly(0,1,2))(tuple(newP[1]))
            newQ = copy(Q)
            newQ.scale_by([1/Q[0][l],1])
            if PK[1][i].abs()<= PK[1][k].abs():
                if k == 0:
                    A = Rx(W.Gpoly(1,0))(tuple(newQ[0]))*PK[1][i]/PK[1][k]
                elif k == 1:
                    A = Rx(W.Gpoly(1,1))(tuple(newQ[0]))*PK[1][i]/PK[1][k]
                else:
                    A = Rx(W.Gpoly(1,2))(tuple(newQ[0]))*PK[1][i]/PK[1][k]
            else:
                if i == 0:
                    A = -Rx(W.Gpoly(1,0))(tuple(newQ[0]))*PK[1][k]/PK[1][i]
                    if k == 1:
                        A = A-Rx(W.Hpoly(1,0,1))(tuple(newQ[0]))
                    else:
                        A = A-Rx(W.Hpoly(1,0,2))(tuple(newQ[0]))
                elif i == 1:
                    A = -Rx(W.Gpoly(1,1))(tuple(newQ[0]))*PK[1][k]/PK[1][i]
                    if k == 0:
                        A = A-Rx(W.Hpoly(1,0,1))(tuple(newQ[0]))
                    else:
                        A = A-Rx(W.Hpoly(1,1,2))(tuple(newQ[0]))
                else:
                    A = -Rx(W.Gpoly(1,2))(tuple(newQ[0]))*PK[1][k]/PK[1][i]
                    if k == 0:
                        A = A-Rx(W.Hpoly(1,0,2))(tuple(newQ[0]))
                    else:
                        A = A-Rx(W.Hpoly(1,1,2))(tuple(newQ[0]))
            #print "A,B:",A,B
            localHeight += beta**(-2*R(e)-1)*R(A.abs()).log() + beta**(-2*R(e))*R(B.abs()).log()
            i = k
            j = l
            newQ.scale_by([1,1/Q[1][k]])
            PK = newQ
        return localHeight

    def canonical_height_plus(self,P,N,badprimes = None,prec = 100):
        r"""
        Evaluates the `+` canonical height function of Call-Silverman
        for ``P`` with ``N`` terms of the series of the local heights.
              Use ``v = 0`` for the archimedean place. Must be over `\ZZ` or `\QQ`.
              ALGORITHM:
              See ``Computing the Canonical Height on K3 Surfaces``, G. Call and J. Silverman,
        Mathematics of Computation, 65(213):259-290, 1996.
              INPUT:
              - ``P`` - a projective point
              - ``N`` - positive integer. number of terms of the series to use
              - ``badprimes` - list of integer primes (where the surface is degenerate) (optional)
              - ``prec`` - float point or p-adic precision, default: 100
              OUTPUT:
              - a real number
              EXAMPLES:::
                  sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ,6)
            sage: L =  (-y0 - y1)*x0 + (-y0*x1 - y2*x2)
            sage: Q = (-y2*y0 - y1^2)*x0^2 + ((-y0^2 - y2*y0 + (-y2*y1 - y2^2))*x1 + (-y0^2 - y2*y1)*x2)*x0 + ((-y0^2 - y2*y0 - y2^2)*x1^2 + (-y2*y0 - y1^2)*x2*x1 + (-y0^2 + (-y1 - y2)*y0)*x2^2)
            sage: X = WehlerK3Surface([L,Q])
            sage: P = X([1,0,-1,1,-1,0]) #order 16
            sage: X.canonical_height_plus(P,5)  # long time
            0.00000000000000000000000000000
                  ::
                  Call-Silverman example
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: P = X(0,1,0,0,0,1)
            sage: X.canonical_height_plus(P,4) # long time
            0.14752753298983071394400412161
        """
        if badprimes == None:
            badprimes = self.degenerate_primes()
        m = 2
        while P[0][m] == 0:
            m = m - 1
            n = 2
        while P[1][n] == 0:
            n = n-1
        h = self.lambda_plus(P,0,N,m,n,prec)
        for p in badprimes:
            h += self.lambda_plus(P,p,N,m,n,prec)
        return(h)

    def canonical_height_minus(self,P,N,badprimes = None,prec = 100):
        r"""
        Evaluates the `+` canonical height function of Call-Silverman
        for ``P`` with ``N`` terms of the series of the local heights.
              Use ``v = 0`` for the archimedean place. Must be over `\ZZ` or `\QQ`.
        ALGORITHM:
            See ``Computing the Canonical Height on K3 Surfaces``, G. Call and J. Silverman,
            Mathematics of Computation, 65(213):259-290, 1996.
        INPUT:
            - ``P`` - a projective point
            - ``N`` - positive integer. number of terms of the series to use
            - ``badprimes` - list of integer primes (where the surface is degenerate) (optional)
            - ``prec`` - float point or p-adic precision, default: 100
        OUTPUT:
            - a real number
        EXAMPLES:::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ,6)
            sage: L =  (-y0 - y1)*x0 + (-y0*x1 - y2*x2)
            sage: Q = (-y2*y0 - y1^2)*x0^2 + ((-y0^2 - y2*y0 + (-y2*y1 - y2^2))*x1 + (-y0^2 - y2*y1)*x2)*x0 + ((-y0^2 - y2*y0 - y2^2)*x1^2 + (-y2*y0 - y1^2)*x2*x1 + (-y0^2 + (-y1 - y2)*y0)*x2^2)
            sage: X = WehlerK3Surface([L,Q])
            sage: P = X([1,0,-1,1,-1,0]) #order 16
            sage: X.canonical_height_minus(P,5)  # long time
            0.00000000000000000000000000000
                
            Call-Silverman example::

            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: P = X(0,1,0,0,0,1)
            sage: X.canonical_height_minus(P,4) # long time
            0.55073705369676788175590206734
        """
        if badprimes == None:
            badprimes = self.degenerate_primes()
        m = 2
        while P[0][m] == 0:
            m = m - 1
            n = 2
        while P[1][n] == 0:
            n = n-1
        h = self.lambda_minus(P,0,N,m,n,prec)
        for p in badprimes:
            h += self.lambda_minus(P,p,N,m,n,prec)
        return(h)

    def canonical_height(self,P,N,badprimes = None,prec = 100):
        r"""
        Evaluates the canonical height for ``P`` with ``N`` terms of the series of the local heights.
        ALGORITHM:
        See ``Computing the Canonical Height on K3 Surfaces``, G. Call and J. Silverman,
        Mathematics of Computation, 65(213):259-290, 1996.
        INPUT:
            - ``P`` - a projective point
            - ``N`` - positive integer. number of terms of the series to use
            - ``badprimes` - list of integer primes (where the surface is degenerate) (optional)
            - ``prec`` - float point or p-adic precision, default: 100
        OUTPUT:
            - a real number
        EXAMPLES::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(QQ,6)
            sage: L =  (-y0 - y1)*x0 + (-y0*x1 - y2*x2)
            sage: Q = (-y2*y0 - y1^2)*x0^2 + ((-y0^2 - y2*y0 + (-y2*y1 - y2^2))*x1 + (-y0^2 - y2*y1)*x2)*x0 + ((-y0^2 - y2*y0 - y2^2)*x1^2 + (-y2*y0 - y1^2)*x2*x1 + (-y0^2 + (-y1 - y2)*y0)*x2^2)
            sage: X = WehlerK3Surface([L,Q])
            sage: P = X([1,0,-1,1,-1,0]) #order 16
            sage: X.canonical_height(P,5)  # long time
            0.00000000000000000000000000000
                  ::
                  Call-Silverman example
            sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: P = X(0,1,0,0,0,1)
            sage: X.canonical_height(P,4)
            0.69826458668659859569990618895
        """
        if badprimes == None:
            badprimes = self.degenerate_primes()
        return(self.canonical_height_plus(P, N,badprimes,precision) + self.canonical_height_minus(P, N,badprimes,precision))

    def Fibers(self,p):
        """
        Returns the fibers of a point on a K3 Surface, will work for nondegenerate fibers only
        INPUT:
           - ``p`` - A list of three integers 
        OUTPUT:
            The corresponding fiber (as a list)
        EXAMPLES:::
            sage: R.<x0,x1,x2,y0,y1,y2> = PolynomialRing(ZZ,6)
            sage: Y = x0*y0+x1*y1-x2*y2
            sage: Z = y0^2*x0*x1 + y0^2*x2^2 - y0*y1*x1*x2 + y1^2*x2*x1 +y2^2*x2^2 + y2^2*x1^2 +y1^2*x2^2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.Fibers([1,0,0])
            TypeError: Fiber is degenerate

            ::

            sage: P.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
            sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
            sage: Y = x0*y0+x1*y1+x2*y2
            sage: X = WehlerK3Surface([Z,Y])
            sage: X.Fibers([0,0,1])
            [(0 : 0 : 1 , 1 : 0 : 0), (0 : 0 : 1 , 0 : 1 : 0)]
        """
        R = self.base_ring()
        Zero = R(0)
        One = R(1)
        P = []
        for i in p:
            j = R(i)
            P.append(j)
        P0 = P+[Zero,Zero,Zero]
        Points = []
        if (self.Gpoly(1,0)(P0)!=0):
            T0 = (self.Hpoly(1,0,1)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,1)(P0))
            T1 = (self.Hpoly(1,0,2)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,2)(P0))
            if (T0.is_square() and T1.is_square()):
                T0 = (self.Hpoly(1,0,1)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,1)(P0)).sqrt()
                T1 = (self.Hpoly(1,0,2)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,2)(P0)).sqrt()
                B1 = (-self.Hpoly(1,0,1)(P0)+T0)/(2*self.Gpoly(1,0)(P0))
                B2 = (-self.Hpoly(1,0,1)(P0)-T0)/(2*self.Gpoly(1,0)(P0))
                C1 = (-self.Hpoly(1,0,2)(P0)+T1)/(2*self.Gpoly(1,0)(P0))
                C2 = (-self.Hpoly(1,0,2)(P0)-T1)/(2*self.Gpoly(1,0)(P0))
                Points.append(P+[One,B1,C1])
                Points.append(P+[One,B2,C1])
                Points.append(P+[One,B1,C2])
                Points.append(P+[One,B2,C2])
            else:
                return []
            
        elif (self.Gpoly(1,1)(P0)!=0):
            T0 = (self.Hpoly(1,0,1)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,1)(P0))
            T1 = (self.Hpoly(1,1,2)(P0)**2-4*self.Gpoly(1,1)(P0)*self.Gpoly(1,2)(P0))
            if (T0.is_square() and T1.is_square()):
                T0 = (self.Hpoly(1,0,1)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,1)(P0)).sqrt()
                T1 = (self.Hpoly(1,1,2)(P0)**2-4*self.Gpoly(1,1)(P0)*self.Gpoly(1,2)(P0)).sqrt()
                A1 = (-self.Hpoly(1,0,1)(P0)+T0)/(2*self.Gpoly(1,1)(P0))
                A2 = (-self.Hpoly(1,0,1)(P0)-T0)/(2*self.Gpoly(1,1)(P0))    
                C1 = (-self.Hpoly(1,1,2)(P0)+T1)/(2*self.Gpoly(1,1)(P0))
                C2 = (-self.Hpoly(1,1,2)(P0)-T1)/(2*self.Gpoly(1,1)(P0))
                Points.append(P+[A1,One,C1])
                Points.append(P+[A1,One,C2])
                Points.append(P+[A2,One,C1])
                Points.append(P+[A2,One,C2])
            else:
                return []
        elif(self.Gpoly(1,2)(P0)!=0):
            T0 = (self.Hpoly(1,0,2)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,2)(P0))
            T1 = (self.Hpoly(1,1,2)(P0)**2-4*self.Gpoly(1,1)(P0)*self.Gpoly(1,2)(P0))
            if (T0.is_square() and T1.is_square()):
                T0 = (self.Hpoly(1,0,2)(P0)**2-4*self.Gpoly(1,0)(P0)*self.Gpoly(1,2)(P0)).sqrt()
                T1 = (self.Hpoly(1,1,2)(P0)**2-4*self.Gpoly(1,1)(P0)*self.Gpoly(1,2)(P0)).sqrt()
                A1 = (-self.Hpoly(1,0,2)(P0)+T0)/(2*self.Gpoly(1,2)(P0))
                A2 = (-self.Hpoly(1,0,2)(P0)-T0)/(2*self.Gpoly(1,2)(P0))    
                B1 = (-self.Hpoly(1,1,2)(P0)+T1)/(2*self.Gpoly(1,2)(P0))
                B2 = (-self.Hpoly(1,1,2)(P0)-T1)/(2*self.Gpoly(1,2)(P0))
                Points.append(P+[A1,B1,One])
                Points.append(P+[A1,B2,One])
                Points.append(P+[A2,B1,One])
                Points.append(P+[A2,B2,One])
            else:
                return []
        elif(self.Hpoly(1,0,1)(P0) != 0):
            Points.append(P+[Zero,One,Zero])
            Points.append(P+[-self.Hpoly(1,0,1)(P0),Zero,-self.Hpoly(1,1,2)(P0)])
            Points.append(P+[One,Zero,Zero])
            Points.append(P+[Zero,-self.Hpoly(1,0,1)(P0),-self.Hpoly(1,0,2)(P0)])
        elif(self.Hpoly(1,0,2)(P0)!=0):
            Points.append(P+[Zero,Zero,One])
            Points.append(P+[-self.Hpoly(1,0,2)(P0),-self.Hpoly(1,1,2)(P0),Zero])
        elif(self.Hpoly(1,1,2)(P0) != 0):
            #make check = False after verification
            Points.append(P+[Zero,Zero,One])
            Points.append(P+[Zero,One,Zero])
        else:
            raise TypeError, "Fiber is degenerate"

        Fibers = []
        for x in Points:
            if (self.L(x) == 0) and (self.Q(x) == 0):
                Y = self.point(x,False)
                Fibers.append(Y)
        AS = self.ambient_space()
        Alist = [AS(list(P)) for P in Fibers]
        for c in Alist:
            l = Alist.count(c)
            for i in range(l-1):
                Alist.remove(c)
        Fibers = [self(list(point)) for point in Alist]
        return Fibers
    def orbit_phi(self,P,N):
        if (isinstance(N,(list,tuple)) == False):
            N = [0,N]
        try:
            N[0] = ZZ(N[0])
            N[1] = ZZ(N[1])
        except TypeError:
            raise TypeError, "Orbit bounds must be integers"
        if N[0]<0 or N[1] <0:
            raise TypeError, "Orbit bounds must be non-negative"
        if N[0] > N[1]:
            return([])
        Q = copy(P)
        Q.normalize_coordinates()
        for i in range(1,N[0]+1):
            Q = self.phi(Q)
            Q.normalize_coordinates()
        Orb = [Q]
        for i in range(N[0]+1,N[1]+1):
            Q = self.phi(Q)
            Q.normalize_coordinates()
            Orb.append(Q)
        return(Orb)

    def orbit_psi(self,P,N):
        if (isinstance(N,(list,tuple)) == False):
            N = [0,N]
        try:
            N[0] = ZZ(N[0])
            N[1] = ZZ(N[1])
        except TypeError:
            raise TypeError, "Orbit bounds must be integers"
        if N[0]<0 or N[1] <0:
            raise TypeError, "Orbit bounds must be non-negative"
        if N[0] > N[1]:
            return([])
        Q = copy(P)
        Q.normalize_coordinates()
        for i in range(1,N[0]+1):
            Q = self.psi(Q)
            Q.normalize_coordinates()
        Orb = [Q]
        for i in range(N[0]+1,N[1]+1):
            Q = self.psi(Q)
            Q.normalize_coordinates()
            Orb.append(Q)
        return(Orb)

    def is_isomorphic(self,other):
        """
              EXAMPLES::
                sage: PP.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2,2],QQ)
                sage: Z = x0^2*y0^2+3*x0*x1*y0^2+x1^2*y0^2+4*x0^2*y0*y1+3*x0*x1*y0*y1-2*x2^2*y0*y1-x0^2*y1^2+2*x1^2*y1^2-x0*x2*y1^2-4*x1*x2*y1^2+5*x0*x2*y0*y2-4*x1*x2*y0*y2+7*x0^2*y1*y2+4*x1^2*y1*y2+x0*x1*y2^2+3*x2^2*y2^2
                sage: Y = x0*y0+x1*y1+x2*y2
                sage: X = WehlerK3Surface([Z,Y])
                sage: W = WehlerK3Surface([Z+Y^2,Y])
                sage: X.is_isomorphic(W)
                True
        """
        return self.defining_ideal() == other.defining_ideal()
    def is_symmetric_orbit(self,orbit):
        sym = False
        i = 0
        while i < len(orbit)-1 and sym == False: 
            P = orbit[i]
            Q = orbit[i+1]
            if P[0] == Q[0] or P[1] == Q[1]:
                sym = True
            i += 1
        return(sym)

class WehlerK3Surface_field( WehlerK3Surface_ring):
    pass

class WehlerK3Surface_finite_field( WehlerK3Surface_field):
    def FiberCounts( self):         
        def getPx1():
            return ([x, y, 1] for x in self.base_ring() for y in self.base_ring())
        def getPx2():
            return ([x, 1, 0] for x in self.base_ring())
        Count = 0
        Xpoint = [1,0,0]
        Ypoint = [1,0,0]
                   #Create all possible Px1 Values
        for i in getPx1():
            for j in getPx1():
                A = i+j
                if(self.L(A) == 0 and self.Q(A) == 0):
                    Count+= 1
            for k in getPx2():
                A = i+k
                if(self.L(A) == 0 and self.Q(A) == 0):
                    Count+= 1
            B = i+Ypoint
            if(self.L(B) == 0 and self.Q(B) == 0):
                Count += 1
        #Create all possible Px2 Values
        for i in getPx2():
            for j in getPx1():
                A = i+j
                if (self.L(A) == 0 and self.Q(A) == 0):
                    Count += 1
            for k in getPx2():
                A = i+k
                if (self.L(A) == 0 and self.Q(A) == 0):
                    Count += 1
            B = i+Ypoint
            if (self.L(B) == 0 and self.Q(B) == 0):
                Count += 1
        #Create all Xpoint values
        for j in getPx1():
            A = Xpoint+j
            if (self.L(A) == 0 and self.Q(A) == 0):
                Count += 1
        for k in getPx2():
            B = Xpoint+k
            if (self.L(B) == 0 and self.Q(B) == 0):
                Count += 1
            C = Xpoint+Ypoint
        if (self.L(C) == 0 and self.Q(C) == 0):
            Count+= 1
        return Count
    def Surgraph(self):
        #need __cmp__ for non str
        #does not work for singular surfaces
        def getPx1():
            return ([x, y, 1] for x in self.base_ring() for y in self.base_ring())
        def getPx2():
            return ([x, 1, 0] for x in self.base_ring())
        R = self.base_ring()
        Xpoint = [R(1),R(0),R(0)]
        V = []
        E = []
        p = self.base_ring().characteristic()
        for P in getPx1():
            try:
                T = self.Fibers(P)
                for t in T:
                    t[1].normalize_coordinates()
                    V.append(t)
                    Q = self.phi(t)
                    Q.normalize_coordinates()
                    E.append([Q])
            except TypeError:  #degenerate fiber
                for PP in getPx1():
                    A = P+PP
                    if(self.L(A) == 0 and self.Q(A) == 0):
                        A = self(A)
                        V.append(A)
                        Q = self.phi(A)
                        Q.normalize_coordinates()
                        E.append([Q])
                for PP in getPx2():
                    A = P+PP
                    if(self.L(A) == 0 and self.Q(A) == 0):
                        A = self(A)
                        V.append(A)
                        Q = self.phi(A)
                        Q.normalize_coordinates()
                        E.append([Q])
                A = P+Xpoint
                if(self.L(A) == 0 and self.Q(A) == 0):
                    A = self(A)
                    V.append(A)
                    Q = self.phi(A)
                    Q.normalize_coordinates()
                    E.append([Q])
        for P in getPx2():
            try:
                T = self.Fibers(P)
                for t in T:
                    t[1].normalize_coordinates()
                    V.append(t)
                    Q = self.phi(t)
                    Q.normalize_coordinates()
                    E.append([Q])
            except TypeError:  #degenerate fiber
                for PP in getPx1():
                    A = P+PP
                    if(self.L(A) == 0 and self.Q(A) == 0):
                        A = self(A)
                        V.append(A)
                        Q = self.phi(A)
                        Q.normalize_coordinates()
                        E.append([Q])
                for PP in getPx2():
                    A = P+PP
                    if(self.L(A) == 0 and self.Q(A) == 0):
                        A = self(A)
                        V.append(A)
                        Q = self.phi(A)
                        Q.normalize_coordinates()
                        E.append([Q])
                A = P+Xpoint
                if(self.L(A) == 0 and self.Q(A) == 0):
                    A = self(A)
                    V.append(A)
                    Q = self.phi(A)
                    Q.normalize_coordinates()
                    E.append([Q])
        try:
            T = self.Fibers(Xpoint)
            for t in T:
                t[1].normalize_coordinates()
                V.append(t)
                Q = self.phi(t)
                Q.normalize_coordinates()
                E.append([Q])
        except TypeError:  #degenerate fiber
            for Q in getPx1():
                A = Xpoint+Q
                if(self.L(A) == 0 and self.Q(A) == 0):
                    A = self(A)
                    V.append(A)
                    Q = self.phi(A)
                    Q.normalize_coordinates()
                    E.append([Q])
            for Q in getPx2():
                A = Xpoint+Q
                if(self.L(A) == 0 and self.Q(A) == 0):
                    A = self(A)
                    V.append(A)
                    Q = self.phi(A)
                    Q.normalize_coordinates()
                    E.append([Q])
            A = Xpoint+Xpoint
            if(self.L(A) == 0 and self.Q(A) == 0):
                A = self(A)
                V.append(A)
                Q = self.phi(A)
                Q.normalize_coordinates()
                E.append([Q])
        from sage.graphs.digraph import DiGraph
        DiGraph(dict(zip(V,E)), loops = True)
        return(g)

    def CycleCount(self,Data,x):
        Count = 0
        Counts = []
        for p,L in Data:
            Total = 0
            for j in L:
                Total = Total+j
                if j <=  (p*x):
                    Count += j
            Prob = Count/Total
            Counts.append(Prob)
        return Counts