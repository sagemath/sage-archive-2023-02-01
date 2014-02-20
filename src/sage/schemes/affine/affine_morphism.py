r"""
Morphisms on affine varieties

A morphism of schemes determined by rational functions that define
what the morphism does on points in the ambient affine space.


AUTHORS:

- David Kohel, William Stein

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (2013-03) iteration functionality and new directory structure
  for affine/projective
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
from sage.misc.misc                import prod
from sage.rings.all                import Integer, moebius
from sage.rings.arith              import lcm
from sage.rings.complex_field      import ComplexField
from sage.rings.integer_ring       import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.quotient_ring      import QuotientRing_generic
from sage.rings.real_mpfr          import RealField
from sage.schemes.generic.morphism import SchemeMorphism_polynomial





class SchemeMorphism_polynomial_affine_space(SchemeMorphism_polynomial):
    """
    A morphism of schemes determined by rational functions that define
    what the morphism does on points in the ambient affine space.

    EXAMPLES::

        sage: RA.<x,y> = QQ[]
        sage: A2 = AffineSpace(RA)
        sage: RP.<u,v,w> = QQ[]
        sage: P2 = ProjectiveSpace(RP)
        sage: H = A2.Hom(P2)
        sage: f = H([x, y, 1])
        sage: f
        Scheme morphism:
          From: Affine Space of dimension 2 over Rational Field
          To:   Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (x : y : 1)
    """
    def __init__(self, parent, polys, check=True):
        r"""
        The Python constructor.

        See :class:`SchemeMorphism_polynomial` for details.

        INPUT:

        - ``parent`` -- Hom

        - ``polys`` -- list or tuple of polynomial or rational functions

        - ``check`` -- Boolean

        OUTPUT:

        - :class:`SchemeMorphism_polynomial_affine_space`

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: H=Hom(A,A)
            sage: H([3/5*x^2,y^2/(2*x^2)])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[3/5*x^2, y^2/(2*x^2)]) must be rational functions in
            Multivariate Polynomial Ring in x, y over Integer Ring

        ::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: H=Hom(A,A)
            sage: H([3*x^2/(5*y),y^2/(2*x^2)])
            Scheme endomorphism of Affine Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x, y) to
                    (3*x^2/(5*y), y^2/(2*x^2))


            sage: A.<x,y>=AffineSpace(QQ,2)
            sage: H=Hom(A,A)
            sage: H([3/2*x^2,y^2])
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (3/2*x^2, y^2)


            sage: A.<x,y>=AffineSpace(QQ,2)
            sage: X=A.subscheme([x-y^2])
            sage: H=Hom(X,X)
            sage: H([9/4*x^2,3/2*y])
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2
            over Rational Field defined by:
              -y^2 + x
              Defn: Defined on coordinates by sending (x, y) to
                    (9/4*x^2, 3/2*y)

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: H=Hom(P,P)
            sage: f=H([5*x^3 + 3*x*y^2-y^3,3*z^3 + y*x^2, x^3-z^3])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x0, x1) to
                    ((5*x0^3 + 3*x0*x1^2 - x1^3)/(x0^3 - 1), (x0^2*x1 + 3)/(x0^3 - 1))
        """
        if check:
            if not isinstance(polys, (list, tuple)):
                raise TypeError("polys (=%s) must be a list or tuple"%polys)
            source_ring =parent.domain().ambient_space().coordinate_ring()
            target = parent.codomain().ambient_space()
            if len(polys) != target.ngens():
                raise ValueError("there must be %s polynomials"%target.ngens())
            try:
                polys = [source_ring(poly) for poly in polys]
            except TypeError:
                if all(p.base_ring()==source_ring.base_ring() for p in polys)==False:
                    raise TypeError("polys (=%s) must be rational functions in %s"%(polys,source_ring))
                try:
                    polys = [source_ring(poly.numerator())/source_ring(poly.denominator()) for poly in polys]
                except TypeError:
                    raise TypeError("polys (=%s) must be rational functions in %s"%(polys,source_ring))
            if isinstance(source_ring, QuotientRing_generic):
                polys = [f.lift() for f in polys]
        SchemeMorphism_polynomial.__init__(self, parent,polys, False)

    def homogenize(self,n,newvar='h'):
        r"""
        Return the homogenization of ``self``. If ``self.domain()`` is a subscheme, the domain of
        the homogenized map is the projective embedding of ``self.domain()``

        INPUT:

        - ``newvar`` -- the name of the homogenization variable (only used when ``self.domain()`` is affine space)

        - ``n`` -- the n-th projective embedding into projective space

        OUTPUT:

        - :class:`SchemMorphism_polynomial_projective_space`

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: H=Hom(A,A)
            sage: f=H([(x^2-2)/x^5,y^2])
            sage: f.homogenize(2,'z')
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2*z^5 - 2*z^7 : x^5*y^2 : x^5*z^2)

        ::

            sage: A.<x,y>=AffineSpace(CC,2)
            sage: H=Hom(A,A)
            sage: f=H([(x^2-2)/(x*y),y^2-x])
            sage: f.homogenize(0,'z')
            Scheme endomorphism of Projective Space of dimension 2 over Complex
            Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x*y*z^2 : x^2*z^2 + (-2.00000000000000)*z^4 : x*y^3 - x^2*y*z)

        ::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: X=A.subscheme([x-y^2])
            sage: H=Hom(X,X)
            sage: f=H([9*y^2,3*y])
            sage: f.homogenize(2)
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 2 over Integer Ring defined by:
              -x1^2 + x0*x2
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (9*x0*x2 : 3*x1*x2 : x2^2)

        ::

            sage: R.<t>=PolynomialRing(ZZ)
            sage: A.<x,y>=AffineSpace(R,2)
            sage: H=Hom(A,A)
            sage: f=H([(x^2-2)/y,y^2-x])
            sage: f.homogenize(0,'z')
            Scheme endomorphism of Projective Space of dimension 2 over Univariate
            Polynomial Ring in t over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (y*z^2 : x^2*z + (-2)*z^3 : y^3 - x*y*z)
        """
        A=self.domain()
        B=self.codomain()
        N=A.ambient_space().dimension_relative()
        NB=B.ambient_space().dimension_relative()
        Vars=list(A.ambient_space().variable_names())+[newvar]
        S=PolynomialRing(A.base_ring(),Vars)
        try:
            l=lcm([self[i].denominator() for i in range(N)])
        except Exception:  #no lcm
            l=prod([self[i].denominator() for i in range(N)])

        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_generic import MPolynomialRing_generic
        if self.domain().base_ring()==RealField() or self.domain().base_ring()==ComplexField():
            F=[S(((self[i]*l).numerator())._maxima_().divide(self[i].denominator())[0].sage()) for i in range(N)]
        elif isinstance(self.domain().base_ring(),(PolynomialRing_general,MPolynomialRing_generic)):
            F=[S(((self[i]*l).numerator())._maxima_().divide(self[i].denominator())[0].sage()) for i in range(N)]
        else:
            F=[S(self[i]*l) for i in range(N)]
        F.insert(n,S(l))
        d=max([F[i].degree() for i in range(N+1)])
        F=[F[i].homogenize(newvar)*S.gen(N)**(d-F[i].degree()) for i in range(N+1)]
        from sage.schemes.affine.affine_space import is_AffineSpace
        if is_AffineSpace(A)==True:
            from sage.schemes.projective.projective_space import ProjectiveSpace
            X=ProjectiveSpace(A.base_ring(),NB,Vars)
        else:
            X=A.projective_embedding(n).codomain()
            phi=S.hom(X.ambient_space().gens(),X.ambient_space().coordinate_ring())
            F=[phi(f) for f in F]
        H=Hom(X,X)
        return(H(F))

    def dynatomic_polynomial(self,period):
        r"""
        For a map `f:\mathbb{A}^1 \to \mathbb{A}^1` this function computes the (affine) dynatomic polynomial.
        The dynatomic polynomial is the analog of the cyclotomic polynomial and its roots are the points
        of formal period `n`.

        ALGORITHM:

        Homogenize to a map `f:\mathbb{P}^1 \to \mathbb{P}^1` and compute the dynatomic polynomial there.
        Then, dehomogenize.

        INPUT:

        - ``period`` -- a positive integer or a list/tuple `[m,n]` where `m` is the preperiod and `n` is the period

        OUTPUT:

        - If possible, a single variable polynomial in the coordinate ring of ``self``.
          Otherwise a fraction field element of the coordinate ring of ``self``

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(QQ,2)
            sage: H=Hom(A,A)
            sage: f=H([x^2+y^2,y^2])
            sage: f.dynatomic_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: Does not make sense in dimension >1

            ::

            sage: A.<x>=AffineSpace(ZZ,1)
            sage: H=Hom(A,A)
            sage: f=H([(x^2+1)/x])
            sage: f.dynatomic_polynomial(4)
            2*x^12 + 18*x^10 + 57*x^8 + 79*x^6 + 48*x^4 + 12*x^2 + 1

            ::

            sage: A.<x>=AffineSpace(CC,1)
            sage: H=Hom(A,A)
            sage: f=H([(x^2+1)/(3*x)])
            sage: f.dynatomic_polynomial(3)
            13.0000000000000*x^6 + 117.000000000000*x^4 + 78.0000000000000*x^2 +
            1.00000000000000

            ::

            sage: A.<x>=AffineSpace(QQ,1)
            sage: H=Hom(A,A)
            sage: f=H([x^2-10/9])
            sage: f.dynatomic_polynomial([2,1])
            531441*x^4 - 649539*x^2 - 524880
        """
        if self.domain() != self.codomain():
            raise TypeError("Must have same domain and codomain to iterate")
        from sage.schemes.affine.affine_space import is_AffineSpace
        if is_AffineSpace(self.domain())==False:
            raise NotImplementedError("Not implemented for subschemes")
        if self.domain().dimension_relative()>1:
            raise TypeError("Does not make sense in dimension >1")
        F=self.homogenize(1).dynatomic_polynomial(period)
        if F.denominator()==1:
            R=F.parent()
            Vars=list(R.variable_names())
            Vars.pop()
            S=PolynomialRing(R.base_ring(),Vars)
            phi=R.hom([S.gen(0),1],S)
            return(phi(F))
        else:
            R=F.numerator().parent()
            Vars=list(R.variable_names())
            Vars.pop()
            S=PolynomialRing(R.base_ring(),Vars)
            phi=R.hom([S.gen(0),1],S)
            return(phi(F.numerator())/phi(F.denominator()))

    def nth_iterate_map(self,n):
        r"""
        This function returns the nth iterate of ``self``

        ALGORITHM:

        Uses a form of successive squaring to reducing computations.

        .. TODO:: This could be improved.

        INPUT:

        - ``n`` - a positive integer.

        OUTPUT:

        - A map between Affine spaces

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: H=Hom(A,A)
            sage: f=H([(x^2-2)/(2*y),y^2-3*x])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Affine Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x, y) to
                    ((x^4 - 4*x^2 - 8*y^2 + 4)/(8*y^4 - 24*x*y^2), (2*y^5 - 12*x*y^3
            + 18*x^2*y - 3*x^2 + 6)/(2*y))

        ::

            sage: A.<x>=AffineSpace(QQ,1)
            sage: H=Hom(A,A)
            sage: f=H([(3*x^2-2)/(x)])
            sage: f.nth_iterate_map(3)
            Scheme endomorphism of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    ((2187*x^8 - 6174*x^6 + 6300*x^4 - 2744*x^2 + 432)/(81*x^7 -
            168*x^5 + 112*x^3 - 24*x))

        ::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: X=A.subscheme([x-y^2])
            sage: H=Hom(X,X)
            sage: f=H([9*x^2,3*y])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2
            over Integer Ring defined by:
              -y^2 + x
              Defn: Defined on coordinates by sending (x, y) to
                    (729*x^4, 9*y)
        """
        if self.domain() != self.codomain():
            raise TypeError("Domain and Codomain of function not equal")
        N=self.codomain().ambient_space().dimension_relative()
        F=list(self._polys)
        R=F[0].parent()
        Coord_ring=self.codomain().coordinate_ring()
        D=Integer(n).digits(2)
        if isinstance(Coord_ring, QuotientRing_generic):
            PHI=[Coord_ring.gen(i).lift() for i in range(N)]
        else:
            PHI=[Coord_ring.gen(i) for i in range(N)]
        for i in range(len(D)):
            T=[F[j] for j in range(N)]
            for k in range(D[i]):
                PHI=[PHI[j](T) for j in range(N)]
            if i!=len(D)-1: #avoid extra iterate
                F=[R(F[j](T)) for j in range(N)] #'square'
        H=Hom(self.domain(),self.codomain())
        return(H(PHI))

    def nth_iterate(self,P,n):
        r"""
        Returns the point `self^n(P)`

        INPUT:

        - ``P`` -- a point in ``self.domain()``
        - ``n`` -- a positive integer.

        OUTPUT:

        - a point in ``self.codomain()``

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(QQ,2)
            sage: H=Hom(A,A)
            sage: f=H([(x-2*y^2)/x,3*x*y])
            sage: f.nth_iterate(A(9,3),3)
            (-104975/13123, -9566667)

        ::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: X=A.subscheme([x-y^2])
            sage: H=Hom(X,X)
            sage: f=H([9*y^2,3*y])
            sage: f.nth_iterate(X(9,3),4)
            (59049, 243)

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: A.<x,y>=AffineSpace(FractionField(R),2)
            sage: H=Hom(A,A)
            sage: f=H([(x-t*y^2)/x,t*x*y])
            sage: f.nth_iterate(A(1,t),3)
            ((-t^16 + 3*t^13 - 3*t^10 + t^7 + t^5 + t^3 - 1)/(t^5 + t^3 - 1), -t^9 - t^7 + t^4)

        """
        return(P.nth_iterate(self,n))

    def orbit(self,P,n):
        r"""
        Returns the orbit of `P` by ``self``. If `n` is an integer it returns `[P,self(P),\ldots,self^n(P)]`.

        If `n` is a list or tuple `n=[m,k]` it returns `[self^m(P),\ldots,self^k(P)]`

        INPUT:

        - ``P`` -- a point in ``self.domain()``
        - ``n`` -- a non-negative integer or list or tuple of two non-negative integers

        OUTPUT:

        - a list of points in ``self.codomain()``

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(QQ,2)
            sage: H=Hom(A,A)
            sage: f=H([(x-2*y^2)/x,3*x*y])
            sage: f.orbit(A(9,3),3)
            [(9, 3), (-1, 81), (13123, -243), (-104975/13123, -9566667)]

        ::

            sage: A.<x>=AffineSpace(QQ,1)
            sage: H=Hom(A,A)
            sage: f=H([(x-2)/x])
            sage: f.orbit(A(1/2),[1,3])
            [(-3), (5/3), (-1/5)]

        ::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: X=A.subscheme([x-y^2])
            sage: H=Hom(X,X)
            sage: f=H([9*y^2,3*y])
            sage: f.orbit(X(9,3),(0,4))
            [(9, 3), (81, 9), (729, 27), (6561, 81), (59049, 243)]

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: A.<x,y>=AffineSpace(FractionField(R),2)
            sage: H=Hom(A,A)
            sage: f=H([(x-t*y^2)/x,t*x*y])
            sage: f.orbit(A(1,t),3)
            [(1, t), (-t^3 + 1, t^2), ((-t^5 - t^3 + 1)/(-t^3 + 1), -t^6 + t^3),
            ((-t^16 + 3*t^13 - 3*t^10 + t^7 + t^5 + t^3 - 1)/(t^5 + t^3 - 1), -t^9 -
            t^7 + t^4)]

        """
        return(P.orbit(self,n))

    def global_height(self,prec=None):
        r"""
        Returns the maximum of the heights of the coefficients in any of the coordinate functions of ``self``.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: A.<x>=AffineSpace(QQ,1)
            sage: H=Hom(A,A)
            sage: f=H([1/1331*x^2+4000]);
            sage: f.global_height()
            8.29404964010203

        ::

            sage: R.<x>=PolynomialRing(QQ)
            sage: k.<w>=NumberField(x^2+5)
            sage: A.<x,y>=AffineSpace(k,2)
            sage: H=Hom(A,A)
            sage: f=H([13*w*x^2+4*y, 1/w*y^2]);
            sage: f.global_height(prec=100)
            3.3696683136785869233538671082

        .. TODO::

            add heights to integer.pyx and remove special case
        """
        if self.domain().base_ring() == ZZ:
            if prec is None:
                R = RealField()
            else:
                R = RealField(prec)
            H=R(0)
            for i in range(self.domain().ambient_space().dimension_relative()):
                C=self[i].coefficients()
                h=max([c.abs() for c in C])
                H=max(H,R(h).log())
            return(H)
        H=0
        for i in range(self.domain().ambient_space().dimension_relative()):
            C=self[i].coefficients()
            if C==[]: #to deal with the case self[i]=0
                h=0
            else:
                h=max([c.global_height(prec) for c in C])
            H=max(H,h)
        return(H)

class SchemeMorphism_polynomial_affine_space_field(SchemeMorphism_polynomial_affine_space):
    pass

class SchemeMorphism_polynomial_affine_space_finite_field(SchemeMorphism_polynomial_affine_space_field):

    def cyclegraph(self):
        r"""
        returns Digraph of all orbits of self mod `p`. For subschemes, only points on the subscheme whose
        image are also on the subscheme are in the digraph.

        OUTPUT:

        - a digraph

        EXAMPLES::

            sage: P.<x,y>=AffineSpace(GF(5),2)
            sage: H=Hom(P,P)
            sage: f=H([x^2-y,x*y+1])
            sage: f.cyclegraph()
            Looped digraph on 25 vertices

        ::

            sage: P.<x>=AffineSpace(GF(3^3,'t'),1)
            sage: H=Hom(P,P)
            sage: f=H([x^2-1])
            sage: f.cyclegraph()
            Looped digraph on 27 vertices

        ::

            sage: P.<x,y>=AffineSpace(GF(7),2)
            sage: X=P.subscheme(x-y)
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2])
            sage: f.cyclegraph()
            Looped digraph on 7 vertices
        """
        if self.domain() != self.codomain():
            raise NotImplementedError("Domain and Codomain must be equal")
        V=[]
        E=[]
        from sage.schemes.affine.affine_space import is_AffineSpace
        if is_AffineSpace(self.domain())==True:
            for P in self.domain():
                V.append(str(P))
                Q=self(P)
                E.append([str(Q)])
        else:
            X=self.domain()
            for P in X.ambient_space():
                try:
                    XP=X.point(P)
                    V.append(str(XP))
                    Q=self(XP)
                    E.append([str(Q)])
                except TypeError:  # not on the scheme
                    pass
        from sage.graphs.digraph import DiGraph
        g=DiGraph(dict(zip(V,E)), loops=True)
        return g


