r"""
Morphisms on projective varieties

A morphism of schemes determined by rational functions that define
what the morphism does on points in the ambient projective space.


AUTHORS:

- David Kohel, William Stein

- William Stein (2006-02-11): fixed bug where P(0,0,0) was allowed as
  a projective point.

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (2013-03) iteration functionality and new directory structure
  for affine/projective, height functionality

"""

# Historical note: in trac #11599, V.B. renamed
# * _point_morphism_class -> _morphism
# * _homset_class -> _point_homset

#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.homset        import Hom
from sage.libs.pari.gen            import pari
from sage.misc.cachefunc           import cached_method
from sage.modules.free_module_element import vector
from sage.rings.all                import Integer, moebius
from sage.rings.arith              import gcd, lcm
from sage.rings.complex_field      import ComplexField
from sage.rings.finite_rings.constructor import GF
from sage.rings.integer_ring       import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.quotient_ring      import QuotientRing_generic
from sage.rings.rational_field     import QQ
from sage.rings.real_mpfr          import RealField
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from copy import copy


class SchemeMorphism_polynomial_projective_space(SchemeMorphism_polynomial):
    """
    A morphism of schemes determined by rational functions that define
    what the morphism does on points in the ambient projective space.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: H([y,2*x])
        Scheme endomorphism of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (y : 2*x)

    An example of a morphism between projective plane curves (see :trac:`10297`)::

        sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
        sage: f = x^3+y^3+60*z^3
        sage: g = y^2*z-( x^3 - 6400*z^3/3)
        sage: C = Curve(f)
        sage: E = Curve(g)
        sage: xbar,ybar,zbar = C.coordinate_ring().gens()
        sage: H = C.Hom(E)
        sage: H([zbar,xbar-ybar,-(xbar+ybar)/80])
        Scheme morphism:
          From: Projective Curve over Rational Field defined by x^3 + y^3 + 60*z^3
          To:   Projective Curve over Rational Field defined by -x^3 + y^2*z + 6400/3*z^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (z : x - y : -1/80*x - 1/80*y)

    A more complicated example::

        sage: P2.<x,y,z> = ProjectiveSpace(2,QQ)
        sage: P1 = P2.subscheme(x-y)
        sage: H12 = P1.Hom(P2)
        sage: H12([x^2,x*z, z^2])
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x - y
          To:   Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
              (y^2 : y*z : z^2)

    We illustrate some error checking::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: f = H([x-y, x*y])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - y, x*y]) must be of the same degree

        sage: H([x-1, x*y+x])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - 1, x*y + x]) must be homogeneous

        sage: H([exp(x),exp(y)])
        Traceback (most recent call last):
        ...
        TypeError: polys (=[e^x, e^y]) must be elements of
        Multivariate Polynomial Ring in x, y over Rational Field
    """

    def __init__(self, parent, polys, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_polynomial` for details.

        EXAMPLES::

            sage: P1.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = P1.Hom(P1)
            sage: H([y,2*x])
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (y : 2*x)
        """
        SchemeMorphism_polynomial.__init__(self, parent, polys, check)
        if check:
            # morphisms from projective space are always given by
            # homogeneous polynomials of the same degree
            try:
                d = polys[0].degree()
            except AttributeError:
                polys = [f.lift() for f in polys]
            if not all([f.is_homogeneous() for f in polys]):
                raise  ValueError("polys (=%s) must be homogeneous"%polys)
            degs = [f.degree() for f in polys]
            if not all([d == degs[0] for d in degs[1:]]):
                raise ValueError("polys (=%s) must be of the same degree"%polys)

    def __eq__(self, right):
        """
        Tests the equality of two projective spaces.

        INPUT:

        - ``right`` - a map on projective space

        OUTPUT:

        - Boolean - True if ``self`` and ``right`` define the same projective map. False otherwise.

        EXAMPLES::

            sage: P1.<x,y> = ProjectiveSpace(RR,1)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1==P2
            False

            ::

            sage: R.<x,y> = QQ[]
            sage: P1 = ProjectiveSpace(R)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1==P2
            True
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return False
        else:
            n = len(self._polys)
            for i in range(0,n):
                for j in range(i+1,n):
                    if self._polys[i]*right._polys[j] != self._polys[j]*right._polys[i]:
                        return False
        return True

    def __ne__(self, right):
        """
        Tests the inequality of two projective spaces.

        INPUT:

        - ``right`` -- a map on projective space

        OUTPUT:

        - Boolean -- True if ``self`` and ``right`` define different projective maps. False otherwise.

        EXAMPLES::

            sage: P1.<x,y> = ProjectiveSpace(RR,1)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1!=P2
            True

            ::

            sage: R.<x,y> = QQ[]
            sage: P1 = ProjectiveSpace(R)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1!=P2
            False
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return True
        else:
            n=len(self._polys)
            for i in range(0,n):
                for j in range(i+1,n):
                    if self._polys[i]*right._polys[j] != self._polys[j]*right._polys[i]:
                        return True
        return False

    def scale_by(self, t):
        """
        Scales each coordinates by a factor of `t`.

        A ``TypeError`` occurs if the point is not in the coordinate_ring
        of the parent after scaling.

        INPUT:

        - ``t`` -- a ring element

        OUTPUT:

        - None.

        EXAMPLES::

            sage: A.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(A,A)
            sage: f = H([x^3-2*x*y^2,x^2*y])
            sage: f.scale_by(1/x)
            sage: f
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 - 2*y^2 : x*y)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: H = Hom(P,P)
            sage: f = H([3/5*x^2,6*y^2])
            sage: f.scale_by(5/3*t); f
            Scheme endomorphism of Projective Space of dimension 1 over Univariate
            Polynomial Ring in t over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (t*x^2 : 10*t*y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.scale_by(x-y);f
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Finite Field of size 7 defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x*y^2 - y^3 : x*y^2 - y^3 : x*z^2 - y*z^2)
        """
        if t == 0:
            raise ValueError("Cannot scale by 0")
        R=self.domain().coordinate_ring()
        if isinstance(R, QuotientRing_generic):
            phi=R.coerce_map_from(self.domain().ambient_space().coordinate_ring())
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                self._polys[i]=phi(self._polys[i]*t).lift()
        else:
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                self._polys[i]=R(self._polys[i]*t)

    def normalize_coordinates(self):
        """
        Scales by 1/gcd of the coordinate functions. Also, scales to clear any denominators from the coefficients.
        This is done in place.

        OUTPUT:

        - None.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([5/4*x^3,5*x*y^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : 4*y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^3+x*y^2,x*y^2,x*z^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Finite Field of size 7 defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (2*y^2 : y^2 : z^2)

        .. NOTE:: gcd raises an error if the base_ring does not support gcds.
        """
        GCD = gcd(self[0],self[1])
        index=2
        if self[0].lc()>0 or self[1].lc() >0:
            neg=0
        else:
            neg=1
        N=self.codomain().ambient_space().dimension_relative()+1
        while GCD!=1 and index < N:
            if self[index].lc()>0:
                neg=0
            GCD=gcd(GCD,self[index])
            index+=+1

        if GCD != 1:
            R=self.domain().base_ring()
            if neg == 1:
                self.scale_by(R(-1)/GCD)
            else:
                self.scale_by(R(1)/GCD)
        else:
            if neg == 1:
                self.scale_by(-1)

        #clears any denominators from the coefficients
        LCM = lcm([self[i].denominator() for i in range(N)])
        self.scale_by(LCM)

        #scales by 1/gcd of the coefficients.
        GCD = gcd([self[i].content() for i in range(N)])
        if GCD!=1:
            self.scale_by(1/GCD)

    def dynatomic_polynomial(self, period):
        r"""
        For a map `f:\mathbb{P}^1 \to \mathbb{P}^1` this function computes the dynatomic polynomial.

        The dynatomic polynomial is the analog of the cyclotomic
        polynomial and its roots are the points of formal period `n`.

        ALGORITHM:

        For a positive integer `n`, let `[F_n,G_n]` be the coordinates of the `nth` iterate of `f`.
        Then construct

        .. MATH::

            \Phi^{\ast}_n(f)(x,y) = \sum_{d \mid n} (yF_d(x,y) - xG_d(x,y))^{\mu(n/d)}

        where `\mu` is the Moebius function.

        For a pair `[m,n]`, let `f^m = [F_m,G_m]`. Compute

        .. MATH::

            \Phi^{\ast}_{m,n}(f)(x,y) = \Phi^{\ast}_n(f)(F_m,G_m)/\Phi^{\ast}_n(f)(F_{m-1},G_{m-1})

        REFERENCES:

        .. [Hutz] B. Hutz. Efficient determination of rational preperiodic
           points for endomorphisms of projective space.
           :arxiv:`1210.6246`, 2012.

        .. [MoPa] P. Morton and P. Patel. The Galois theory of periodic points
           of polynomial maps. Proc. London Math. Soc., 68 (1994), 225-263.

        INPUT:

        - ``period`` -- a positive integer or a list/tuple `[m,n]` where `m` is the preperiod and `n` is the period

        OUTPUT:

        - If possible, a two variable polynomial in the coordinate ring of ``self``.
          Otherwise a fraction field element of the coordinate ring of ``self``

        .. TODO::

            Do the division when the base ring is p-adic or a function field
            so that the output is a polynomial.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + 2*y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,x*y])
            sage: f.dynatomic_polynomial(4)
            2*x^12 + 18*x^10*y^2 + 57*x^8*y^4 + 79*x^6*y^6 + 48*x^4*y^8 + 12*x^2*y^10 + y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,3*x*y])
            sage: f.dynatomic_polynomial(3)
            13.0000000000000*x^6 + 117.000000000000*x^4*y^2 +
            78.0000000000000*x^2*y^4 + y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-10/9*y^2,y^2])
            sage: f.dynatomic_polynomial([2,1])
            x^4*y^2 - 11/9*x^2*y^4 - 80/81*y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-29/16*y^2,y^2])
            sage: f.dynatomic_polynomial([2,3])
            x^12 - 95/8*x^10*y^2 + 13799/256*x^8*y^4 - 119953/1024*x^6*y^6 +
            8198847/65536*x^4*y^8 - 31492431/524288*x^2*y^10 +
            172692729/16777216*y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.dynatomic_polynomial([1,2])
            x^2 - x*y

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3-y^3,3*x*y^2])
            sage: f.dynatomic_polynomial([0,4])==f.dynatomic_polynomial(4)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,x*y,z^2])
            sage: f.dynatomic_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: Does not make sense in dimension >1

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(5),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.dynatomic_polynomial(2)
            (x^4*y + (2 + O(5^20))*x^2*y^3 - x*y^4 + (2 + O(5^20))*y^5)/(x^2*y -
            x*y^2 + y^3)

        .. TODO:: It would be nice to get this to actually be a polynomial.

        ::

            sage: L.<t> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(L,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + (t + 1)*y^2

        ::

            sage: K.<c> = PolynomialRing(ZZ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+ c*y^2,y^2])
            sage: f.dynatomic_polynomial([1,2])
            x^2 - x*y + (c + 1)*y^2
       """
        if self.domain() != self.codomain():
            raise TypeError("Must have same domain and codomain to iterate")
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is False:
            raise NotImplementedError("Not implemented for subschemes")
        if self.domain().dimension_relative()>1:
            raise TypeError("Does not make sense in dimension >1")

        if (isinstance(period,(list,tuple)) is False):
            period=[0,period]
        try:
            period[0]=ZZ(period[0])
            period[1]=ZZ(period[1])
        except TypeError:
            raise TypeError("Period and preperiod must be integers")
        if period[1]<=0:
                raise AttributeError("Period must be at least 1")

        if period[0]!=0:
            m=period[0]
            fm=self.nth_iterate_map(m)
            fm1=self.nth_iterate_map(m-1)
            n=period[1]
            PHI=1;
            x=self.domain().gen(0)
            y=self.domain().gen(1)
            F=self._polys
            f=F
            for d in range(1,n+1):
                if n%d == 0:
                    PHI=PHI*((y*F[0]-x*F[1])**moebius(n/d))
                if d !=n: #avoid extra iteration
                    F=[f[0](F[0],F[1]),f[1](F[0],F[1])]
            if m!=0:
                PHI=PHI(fm._polys)/PHI(fm1._polys)
        else:
            PHI=1;
            x=self.domain().gen(0)
            y=self.domain().gen(1)
            F=self._polys
            f=F
            for d in range(1,period[1]+1):
                if period[1]%d == 0:
                    PHI=PHI*((y*F[0]-x*F[1])**moebius(period[1]/d))
                if d !=period[1]: #avoid extra iteration
                    F=[f[0](F[0],F[1]),f[1](F[0],F[1])]
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_generic import MPolynomialRing_generic
        if (self.domain().base_ring() == RealField()
            or self.domain().base_ring() == ComplexField()):
            PHI=PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
        elif isinstance(self.domain().base_ring(),(PolynomialRing_general,MPolynomialRing_generic)):
            from sage.rings.padics.generic_nodes import is_pAdicField, is_pAdicRing
            from sage.rings.function_field.function_field import is_FunctionField
            BR=self.domain().base_ring().base_ring()
            if is_pAdicField(BR) or is_pAdicRing(BR) or is_FunctionField(BR):
                raise NotImplementedError("Not implemented")
            PHI=PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
            #do it again to divide out by denominators of coefficients
            PHI=PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
        if PHI.denominator() == 1:
            PHI=self.coordinate_ring()(PHI)
        return(PHI)

    def nth_iterate_map(self, n):
        r"""
        For a map ``self`` this function returns the nth iterate of ``self`` as a
        function on ``self.domain()``

        ALGORITHM:

        Uses a form of successive squaring to reducing computations.


        .. TODO:: This could be improved.

        INPUT:

        - ``n`` -- a positive integer.

        OUTPUT:

        - A map between projective spaces

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 + 2*x^2*y^2 + 2*y^4 : y^4)

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,x*y])
            sage: f.nth_iterate_map(3)
            Scheme endomorphism of Projective Space of dimension 1 over Complex
            Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y) to
                    (x^8 + (-7.00000000000000)*x^6*y^2 + 13.0000000000000*x^4*y^4 +
            (-7.00000000000000)*x^2*y^6 + y^8 : x^7*y + (-4.00000000000000)*x^5*y^3
            + 4.00000000000000*x^3*y^5 - x*y^7)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,x*y,z^2+x^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 - 3*x^2*y^2 + y^4 : x^3*y - x*y^3 : 2*x^4 - 2*x^2*y^2 + y^4
            + 2*x^2*z^2 + z^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x*z-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,x*z,z^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Rational Field defined by:
              -y^2 + x*z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 : x^2*z^2 : z^4)
        """
        if self.domain() != self.codomain():
            raise TypeError("Domain and Codomain of function not equal")
        if n <0:
            raise TypeError("Iterate number must be a positive integer")
        N=self.codomain().ambient_space().dimension_relative()+1
        F=list(self._polys)
        D=Integer(n).digits(2)  #need base 2
        Coord_ring=self.codomain().coordinate_ring()
        if isinstance(Coord_ring, QuotientRing_generic):
            PHI=[Coord_ring.gen(i).lift() for i in range(N)]
        else:
            PHI=[Coord_ring.gen(i) for i in range(N)]
        for i in range(len(D)):
            T=tuple([F[j] for j in range(N)])
            for k in range(D[i]):
                PHI=[PHI[j](T) for j in range(N)]
            if i!=len(D)-1: #avoid extra iterate
                F=[F[j](T) for j in range(N)] #'square'
        H=Hom(self.domain(),self.codomain())
        return(H(PHI))

    def nth_iterate(self, P, n, normalize=False):
        r"""
        For a map ``self`` and a point `P` in ``self.domain()``
        this function returns the nth iterate of `P` by ``self``.

        If ``normalize`` is ``True``, then the coordinates are
        automatically normalized.

        .. TODO:: Is there a more efficient way to do this?

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        - ``n`` -- a positive integer.

        - ``normalize`` - Boolean (optional Default: ``False``)

        OUTPUT:

        - A point in ``self.codomain()``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q,4)
            (32768 : 32768)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q,4,1)
            (1 : 1)

        Is this the right behavior? ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2,2*y^2,z^2-x^2])
            sage: Q = P(2,7,1)
            sage: f.nth_iterate(Q,2)
            (-16/7 : -2744 : 1)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,(2-t)*y^2,z^2])
            sage: Q = P(2+t,7,t)
            sage: f.nth_iterate(Q,2)
            (t^4 + 2507*t^3 - 6787*t^2 + 10028*t + 16 : -2401*t^3 + 14406*t^2 -
            28812*t + 19208 : t^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.nth_iterate(X(2,2,3),3)
            (256 : 256 : 6561)

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3-2*x*y^2 - c*y^3,x*y^2])
            sage: f.nth_iterate(P(c,1),2)
            ((c^6 - 9*c^4 + 25*c^2 - c - 21)/(c^2 - 3) : 1)
        """
        return(P.nth_iterate(self,n,normalize))

    def degree(self):
        r"""
        This function returns the degree of ``self``.

        The degree is defined as the degree of the homogeneous
        polynomials that are the coordinates of ``self``.

        OUTPUT:

        - A positive integer

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.degree()
            2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(CC,2)
            sage: H = Hom(P,P)
            sage: f = H([x^3+y^3,y^2*z,z*x*y])
            sage: f.degree()
            3

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,(2-t)*y^2,z^2])
            sage: f.degree()
            2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.degree()
            2
        """
        return(self._polys[0].degree())

    def dehomogenize(self, n):
        r"""
        Returns the standard dehomogenization at the nth coordinate `(\frac{self[0]}{self[n]},\frac{self[1]}{self[n]},...)`.

        Note that the new function is defined over the fraction field
        of the base ring of ``self``.

        INPUT:

        - ``n`` -- a nonnegative integer

        OUTPUT:

        - :class:`SchemeMorphism_polynomial_affine_space` (on nth affine patch)

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.dehomogenize(0)
            Scheme endomorphism of Affine Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x) to
                    (x^2/(x^2 + 1))

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2-z^2,2*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1/2*x0^2 + 1/2*x1^2, 1/2*x1^2 - 1/2)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,t*y^2-z^2,t*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Fraction Field
            of Univariate Polynomial Ring in t over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1/t*x0^2 + x1^2, x1^2 - 1/t)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,x*z])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2
            over Integer Ring defined by:
              x0^2 - x1^2
              Defn: Defined on coordinates by sending (x0, x1) to
                    (x1^2/x0, x1^2/x0)
        """
        PS=self.domain()
        A=PS.ambient_space()
        if self._polys[n].substitute({A.gen(n):1}) == 0:
            raise ValueError("Can't dehomogenize at 0 coordinate.")
        else:
            Aff=PS.affine_patch(n)
            S=Aff.ambient_space().coordinate_ring()
            R=A.coordinate_ring()
            phi=R.hom([S.gen(j) for j in range(0,n)] + [1] + [S.gen(j) for j in range(n,A.dimension_relative())],S)
            F=[]
            for i in range(0,A.dimension_relative()+1):
                if i !=n:
                    F.append(phi(self._polys[i])/phi(self._polys[n]))
            H=Hom(Aff,Aff)
            return(H(F))

    def orbit(self, P, N, **kwds):
        r"""
        Returns the orbit of `P` by ``self``. If `n` is an integer it returns `[P,self(P),\ldots,self^n(P)]`.
        If `n` is a list or tuple `n=[m,k]` it returns `[self^m(P),\ldots,self^k(P)]`.
        Automatically normalize the points if ``normalize=True``. Perform the checks on point initialize if
        ``check=True``

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        - ``n`` -- a non-negative integer or list or tuple of two non-negative integers

        kwds:

        - ``check`` -- boolean (optional - default: ``True``)

        - ``normalize`` -- boolean (optional - default: ``False``)


        OUTPUT:

        - a list of points in ``self.codomain()``

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2-z^2,2*z^2])
            sage: f.orbit(P(1,2,1),3)
            [(1 : 2 : 1), (5 : 3 : 2), (34 : 5 : 8), (1181 : -39 : 128)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2-z^2,2*z^2])
            sage: f.orbit(P(1,2,1),[2,4])
            [(34 : 5 : 8), (1181 : -39 : 128), (1396282 : -14863 : 32768)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,x*z])
            sage: f.orbit(X(2,2,3),3,normalize=True)
            [(2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.orbit(P.point([1,2],False),4,check=False)
            [(1 : 2), (5 : 4), (41 : 16), (1937 : 256), (3817505 : 65536)]

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+c*y^2,y^2])
            sage: f.orbit(P(0,1),3)
            [(0 : 1), (c : 1), (c^2 + c : 1), (c^4 + 2*c^3 + c^2 + c : 1)]
        """
        return(P.orbit(self,N,**kwds))

    @cached_method
    def is_morphism(self):
        r"""
        returns ``True`` if self is a morphism (no common zero of defining polynomials).

        The map is a morphism if and only if the ideal generated by
        the defining polynomials is the unit ideal.

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.is_morphism()
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(RR,2)
            sage: H = Hom(P,P)
            sage: f = H([x*z-y*z,x^2-y^2,z^2])
            sage: f.is_morphism()
            False

        ::

            sage: R.<t> = PolynomialRing(GF(5))
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x*z-t*y^2,x^2-y^2,t*z^2])
            sage: f.is_morphism()
            True
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is False or is_ProjectiveSpace(self.codomain()) is False:
            raise NotImplementedError
        R=self.coordinate_ring()
        F=self._polys
        if R.base_ring().is_field():
            J=R.ideal(F)
        else:
            S=PolynomialRing(R.base_ring().fraction_field(),R.gens(),R.ngens())
            J=S.ideal([S.coerce(F[i]) for i in range(R.ngens())])
        if J.dimension()>0:
            return False
        else:
            return True

    def resultant(self, normalize=False):
        r"""
        Computes the resultant of the defining polynomials of ``self`` if ``self`` is a map on the projective line.

        If ``normalize`` is ``True``, then first normalize the coordinate
        functions with :meth:`normalize_coordinates`.

        INPUT:

        - ``normalize`` -- Boolean (optional - default: ``False``)

        OUTPUT:

        - an element of ``self.codomain().base_ring()``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,6*y^2])
            sage: f.resultant()
            36

        ::

            sage: R.<t> = PolynomialRing(GF(17))
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: H = Hom(P,P)
            sage: f = H([t*x^2+t*y^2,6*y^2])
            sage: f.resultant()
            2*t^2
        """
        if self.domain().dimension_relative() > 1:
            raise TypeError("Only for dimension 1, use self.primes_of_bad_reduction() to get bad primes")
        if normalize is True:
            F=copy(self)
            F.normalize_coordinates()
        else:
            F=self
        x=self.domain().gen(0)
        y=self.domain().gen(1)
        d=self.degree()
        f=F[0].substitute({y:1})
        g=F[1].substitute({y:1})
        res=(f.lc()**(d-g.degree())*g.lc()**(d-f.degree())*pari(f).polresultant(g,x))
        return(self.codomain().base_ring()(res))

    @cached_method
    def primes_of_bad_reduction(self, check=True):
        r"""
        Determines the primes of bad reduction for a map `self: \mathbb{P}^N \to \mathbb{P}^N`
        defined over `\ZZ` or `\QQ`.

        If ``check`` is ``True``, each prime is verified to be of bad reduction.

        ALGORITHM:

        `p` is a prime of bad reduction if and only if the defining
        polynomials of self have a common zero. Or stated another way,
        `p` is a prime of bad reducion if and only if the radical of
        the ideal defined by the defining polynomials of self is not
        `(x_0,x_1,\ldots,x_N)`.  This happens if and only if some
        power of each `x_i` is not in the ideal defined by the
        defining polynomials of self. This last condition is what is
        checked. The lcm of the coefficients of the monomials `x_i` in
        a groebner basis is computed. This may return extra primes.

        INPUT:

        - ``check`` -- Boolean (optional - default: ``True``)

        OUTPUT:

        - a list of integer primes.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/3*x^2+1/2*y^2,y^2])
            sage: print f.primes_of_bad_reduction()
            [2, 3]

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ,3)
            sage: H = Hom(P,P)
            sage: f = H([12*x*z-7*y^2,31*x^2-y^2,26*z^2,3*w^2-z*w])
            sage: f.primes_of_bad_reduction()
            [2, 3, 7, 13, 31]

        This is an example where check=False returns extra primes::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([3*x*y^2 + 7*y^3 - 4*y^2*z + 5*z^3, -5*x^3 + x^2*y + y^3 + 2*x^2*z, -2*x^2*y + x*y^2 + y^3 - 4*y^2*z + x*z^2])
            sage: f.primes_of_bad_reduction(False)
            [2, 5, 37, 2239, 304432717]
            sage: f.primes_of_bad_reduction()
            [5, 37, 2239, 304432717]
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise TypeError("Must be ZZ or QQ")
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is False or is_ProjectiveSpace(self.codomain()) is False:
            raise NotImplementedError
        R=self.coordinate_ring()
        F=self._polys
        if R.base_ring().is_field():
            J=R.ideal(F)
        else:
            S=PolynomialRing(R.base_ring().fraction_field(),R.gens(),R.ngens())
            J=S.ideal([S.coerce(F[i]) for i in range(R.ngens())])
        if J.dimension()>0:
            raise TypeError("Not a morphism.")
        #normalize to coefficients in the ring not the fraction field.
        F=[F[i]*lcm([F[j].denominator() for j in range(len(F))]) for i in range(len(F))]

        #move the ideal to the ring of integers
        if R.base_ring().is_field():
            S=PolynomialRing(R.base_ring().ring_of_integers(),R.gens(),R.ngens())
            F=[F[i].change_ring(R.base_ring().ring_of_integers()) for i in range(len(F))]
            J=S.ideal(F)
        else:
            J=R.ideal(F)
        GB=J.groebner_basis()
        badprimes=[]

        #get the primes dividing the coefficients of the monomials x_i^k_i
        for i in range(len(GB)):
            LT=GB[i].lt().degrees()
            power=0
            for j in range(R.ngens()):
                if LT[j]!=0:
                    power+=1
            if power == 1:
                badprimes=badprimes+GB[i].lt().coefficients()[0].support()
        badprimes=list(set(badprimes))
        badprimes.sort()

        #check to return only the truly bad primes
        if check == True:
            index=0
            while index < len(badprimes):  #figure out which primes are really bad primes...
                S=PolynomialRing(GF(badprimes[index]),R.gens(),R.ngens())
                J=S.ideal([S.coerce(F[j]) for j in range(R.ngens())])
                if J.dimension() == 0:
                    badprimes.pop(index)
                else:
                    index+=1
        return(badprimes)

    def conjugate(self, M):
        r"""
        Conjugates ``self`` by ``M``, i.e. `M^{-1} \circ f \circ M`.

        If possible the map will be defined over the same space as
        ``self``. Otherwise, will try to coerce to the base_ring of
        ``M``.

        INPUT:

        - ``M`` -- a square invertible matrix

        OUTPUT:

        - a map from ``self.domain()`` to ``self.codomain()``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.conjugate(matrix([[1,2],[0,1]]))
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 4*x*y + 3*y^2 : y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1)
            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3+y^3,y^3])
            sage: f.conjugate(matrix([[i,0],[0,-i]]))
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (-x^3 + y^3 : -y^3)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2,y*z])
            sage: f.conjugate(matrix([[1,2,3],[0,1,2],[0,0,1]]))
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 + 4*x*y + 3*y^2 + 6*x*z + 9*y*z + 7*z^2 : y^2 + 2*y*z : y*z + 2*z^2)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.conjugate(matrix([[2,0],[0,1/2]]))
            Scheme endomorphism of Projective Space of dimension 1 over Multivariate
            Polynomial Ring in x, y over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x^2 + 1/8*y^2 : 1/2*y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/3*x^2+1/2*y^2,y^2])
            sage: f.conjugate(matrix([[i,0],[0,-i]]))
            Scheme endomorphism of Projective Space of dimension 1 over Multivariate
            Polynomial Ring in x, y over Number Field in i with defining polynomial
            x^2 + 1
              Defn: Defined on coordinates by sending (x : y) to
                    ((1/3*i)*x^2 + (1/2*i)*y^2 : (-i)*y^2)
        """

        if M.is_square() == 1 and M.determinant() != 0 and  M.ncols() == self.domain().ambient_space().dimension_relative()+1:
            X=M*vector(self[0].parent().gens())
            F=vector(self._polys)
            F=F(list(X))
            N=M.inverse()
            F=N*F
            R=self.codomain().ambient_space().coordinate_ring()
            try:
                F = [R(f) for f in F]
                PS=self.codomain()
            except TypeError: #no longer defined over same ring
                R=R.change_ring(M.base_ring())
                F = [R(f) for f in F]
                PS=self.codomain().change_ring(R)
            H=Hom(PS,PS)
            return(H(F))
        else:
            raise TypeError("matrix must be invertible and size dimension +1 ")

    def green_function(self, P, v, **kwds):
        r"""
        Evaluates the local Green's function at the place ``v`` for ``P`` with ``N`` terms of the
        series or, in dimension 1, to within a given error bound.

        Use ``v=0`` for the archimedean place. Must be over `\ZZ` or `\QQ`.

        ALGORITHM:

        See Exercise 5.29 and Figure 5.6 of ``The Arithmetic of Dynamics Systems``, Joseph H. Silverman, Springer, GTM 241, 2007.

        INPUT:

        - ``P`` - a projective point

        - ``v`` - non-negative integer. a place, use v=0 for the archimedean place

        kwds:

        - ``N`` - positive integer. number of terms of the series to use

        - ``prec`` - positive integer, float point or p-adic precision, default: 100

        - ``error_bound`` - a positive real number

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,x*y])
            sage: f.green_function(P.point([5,2],False),0,N=30)
            1.7315451844777407992085512000
            sage: f.green_function(P.point([2,1],False),0,N=30)
            0.86577259223181088325226209926
            sage: f.green_function(P.point([1,1],False),0,N=30)
            0.43288629610862338612700146098
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise TypeError("Must be ZZ or QQ")
        return(P.green_function(self,v,**kwds))

    def canonical_height(self, P, **kwds):
        r"""
        Evaluates the canonical height of ``P`` with respect to ``self``. Must be over `\ZZ` or `\QQ`.

        Specify either the number of terms of the series to evaluate
        or, in dimension 1, the error bound required.

        ALGORITHM:

        The sum of the Green's function at the archimedean place and the places of bad reduction.

        INPUT:

        - ``P`` -- a projective point

        kwds:

        - ``badprimes`` - a list of primes of bad reduction

        - ``N`` - positive integer. number of terms of the series to use in the local green functions

        - ``prec`` - positive integer, float point or p-adic precision, default: 100

        - ``error_bound`` - a positive real number

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,2*x*y]);
            sage: f.canonical_height(P.point([5,4]),error_bound=0.001)
            2.1970553519503404898926835324
            sage: f.canonical_height(P.point([2,1]),error_bound=0.001)
            1.0984430632822307984974382955

        Notice that preperiodic points may not be exactly 0::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-29/16*y^2,y^2]);
            sage: f.canonical_height(P.point([1,4]),N=60)
            1.2024186864216154694752186858e-18

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x^2-y^2);
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,4*z^2]);
            sage: Q = X([4,4,1])
            sage: f.canonical_height(Q,badprimes=[2])
            0.0013538030870311431824555314882
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise TypeError("Must be ZZ or QQ")
        return(P.canonical_height(self,**kwds))

    def global_height(self, prec=None):
        r"""
        Returns the maximum of the heights of the coefficients in any of the coordinate functions of ``self``.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/1331*x^2+1/4000*y^2,210*x*y]);
            sage: f.global_height()
            8.29404964010203

        This function does not automatically normalize::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([4*x^2+100*y^2,210*x*y,10000*z^2]);
            sage: f.global_height()
            9.21034037197618
            sage: f.normalize_coordinates()
            sage: f.global_height()
            8.51719319141624

        ::

            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2-2)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O,1)
            sage: H = Hom(P,P)
            sage: f = H([2*x^2 + 3*O(w)*y^2,O(w)*y^2])
            sage: f.global_height()
            1.44518587894808

        .. TODO:: add heights to integer.pyx and remove special case
        """
        if self.domain().base_ring() == ZZ:
            if prec is None:
                R = RealField()
            else:
                R = RealField(prec)
            H=R(0)
            for i in range(self.domain().ambient_space().dimension_relative()+1):
                C=self[i].coefficients()
                h=max([c.abs() for c in C])
                H=max(H,R(h).log())
            return(H)
        H=0
        for i in range(self.domain().ambient_space().dimension_relative()+1):
            C=self[i].coefficients()
            h=max([c.global_height(prec) for c in C])
            H=max(H,h)
        return(H)

class SchemeMorphism_polynomial_projective_space_field(SchemeMorphism_polynomial_projective_space):
    """
    Place holder
    """
    pass

class SchemeMorphism_polynomial_projective_space_finite_field(SchemeMorphism_polynomial_projective_space_field):

    def orbit_structure(self, P):
        r"""
        Every point is preperiodic over a finite field. This funtion returns the pair `[m,n]` where `m` is the
        preperiod and `n` the period of the point ``P`` by ``self``.

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        OUTPUT:

        - a list `[m,n]` of integers

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2,z^2 + y*z])
            sage: f.orbit_structure(P(2,1,2))
            [0, 6]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.orbit_structure(X(1,1,2))
            [0, 2]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.orbit_structure(P(3,4))
            [2, 3]
        """
        return(P.orbit_structure(self))

    def cyclegraph(self):
        r"""
        returns Digraph of all orbits of ``self`` mod `p`.

        For subschemes, only points on the subscheme whose image are
        also on the subscheme are in the digraph.

        OUTPUT:

        - a digraph

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.cyclegraph()
            Looped digraph on 14 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5^2,'t'),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2,z^2+y*z])
            sage: f.cyclegraph()
            Looped digraph on 651 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.cyclegraph()
            Looped digraph on 15 vertices
        """
        if self.domain() != self.codomain():
            raise NotImplementedError("Domain and Codomain must be equal")
        V=[]
        E=[]
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is True:
            for P in self.domain():
                V.append(str(P))
                Q=self(P)
                Q.normalize_coordinates()
                E.append([str(Q)])
        else:
            X=self.domain()
            for P in X.ambient_space():
                try:
                    XP=X.point(P)
                    V.append(str(XP))
                    Q=self(XP)
                    Q.normalize_coordinates()
                    E.append([str(Q)])
                except TypeError:  # not a point on the scheme
                    pass
        from sage.graphs.digraph import DiGraph
        g=DiGraph(dict(zip(V,E)), loops=True)
        return g
