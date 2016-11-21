r"""
Polynomial morphisms for products of projective spaces

This class builds on the projective space class and its point and morphism classes.

EXAMPLES::

    sage: P1xP1.<x,y,u,v> = ProductProjectiveSpaces(QQ, [1, 1])
    sage: H = End(P1xP1)
    sage: H([x^2*u, y^2*v, x*v^2, y*u^2])
    Scheme endomorphism of Product of projective spaces P^1 x P^1 over Rational Field
      Defn: Defined by sending (x : y , u : v) to
            (x^2*u : y^2*v , x*v^2 : y*u^2).
"""
#*****************************************************************************
# Copyright (C) 2014 Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.homset import Hom, End
from sage.rings.quotient_ring import QuotientRing_generic
from sage.schemes.generic.morphism import SchemeMorphism_polynomial


class ProductProjectiveSpaces_morphism_ring(SchemeMorphism_polynomial):
    r"""
    The class of morphisms on products of projective spaces.

    The components are projective space morphisms.

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: H = T.Hom(T)
        sage: H([x^2, y^2, z^2, w^2, u^2])
        Scheme endomorphism of Product of projective spaces P^2 x P^1 over Rational Field
          Defn: Defined by sending (x : y : z , w : u) to
                (x^2 : y^2 : z^2 , w^2 : u^2).
    """

    def __init__(self, parent, polys, check=True):
        r"""
        The Python constructor.

        INPUT:

        - ``parent`` -- Hom-set.

        - ``polys`` -- anything that defines a point in the class.

        - ``check`` -- Boolean. Whether or not to perform input checks.
          (Default:`` True``)

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u, y^2*w, z^2*u, w^2, u^2])
            Scheme endomorphism of Product of projective spaces P^2 x P^1 over Rational Field
              Defn: Defined by sending (x : y : z , w : u) to
                    (x^2*u : y^2*w : z^2*u , w^2 : u^2).

        ::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u, y^2*w, z^2*u, w^2, u*z])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[x^2*u, y^2*w, z^2*u, w^2, z*u]) must be
            multi-homogeneous of the same degrees (by component)

        ::

            sage: R.<s,t> = PolynomialRing(QQ)
            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: P.<u,v,w,s,t,r> = ProductProjectiveSpaces([3, 1], QQ)
            sage: H = Hom(Z,P)
            sage: f = H([a^2,b^2,a^2,a*b,a*x,b*z]); f
            Scheme morphism:
              From: Product of projective spaces P^1 x P^2 over Rational Field
              To:   Product of projective spaces P^3 x P^1 over Rational Field
              Defn: Defined by sending (a : b , x : y : z) to
                    (a^2 : b^2 : a^2 : a*b , a*x : b*z).

        ::

            sage: Z.<a,b,c,x,y,z> = ProductProjectiveSpaces([1, 3], QQ)
            sage: P.<u,v,w,s,t,r> = ProductProjectiveSpaces([2, 2], QQ)
            sage: H = Hom(Z,P)
            sage: f = H([a^2,b^2,c^2,x^2,y^2,z^2])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[a^2, b^2, c^2, x^2, y^2, z^2]) must be
            multi-homogeneous of the same degrees (by component)
        """
        if check:
            #check multi-homogeneous
            #if self is a subscheme, we may need the lift of the polynomials
            try:
                polys[0].exponents()
            except AttributeError:
                polys = [f.lift() for f in polys]

            target = parent.codomain().ambient_space()
            dom = parent.domain().ambient_space()
            from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
            if is_ProductProjectiveSpaces(target):
                splitpolys = target._factors(polys)
                for m in range(len(splitpolys)):
                    d = dom._degree(splitpolys[m][0])
                    if not all(d == dom._degree(f) for f in splitpolys[m]):
                        raise  TypeError("polys (=%s) must be multi-homogeneous of the same degrees (by component)"%polys)
            else:
                #we are mapping into some other kind of space
                target._validate(polys)

        SchemeMorphism_polynomial.__init__(self, parent, polys, check)

    def __getitem__(self, i):
        r"""
        Return the ``i``-th coordinate polynomial.

        INPUT:

        - ``i`` -- integer.

        OUTPUT:

        The (multi)-homomgeneous polynomial that is the ``i``-th coordinate.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u, y^2*w, z^2*u, w^2, u^2])
            sage: F[2]
            z^2*u
        """
        return(self._polys[i])

    def _repr_defn(self):
        r"""
        Return a string representation of this morphism.

        OUTPUT: String.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], QQ)
            sage: H = Hom(P,P)
            sage: f = H([x^2, y^2, z, w])
            sage: f._repr_defn()
            'Defined by sending (x : y , z : w) to \n(x^2 : y^2 , z : w).'
        """
        s  = 'Defined by sending '
        s += self.domain().ambient_space()._repr_generic_point()
        s += ' to \n'
        s += self.codomain().ambient_space()._repr_generic_point(self._polys)
        s += '.'
        return s

    def __call__(self, P, check=True):
        r"""
        Make morphisms of products of projective spaces callable.

        INPUT:

        - ``P`` -- a point in the domain.

        - ``check`` -- Boolean - whether or not to perform the input checks
          on the image point (Default: ``True``).

        OUTPUT: The image point in the codomain.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u, y^2*w, z^2*u, w^2, u^2])
            sage: F(T([2, 1, 3, 0, 1]))
            (4/9 : 0 : 1 , 0 : 1)

        ::

            sage: PP.<x,y,z,u,v,w> = ProductProjectiveSpaces(QQ, [1, 1, 1])
            sage: HP = End(PP)
            sage: f = HP([v*x^2,w*y^2,z^2,u^2,v^2,w^2])
            sage: Q = PP([0,1,1,1,1,1])
            sage: f(Q)
            (0 : 1 , 1 : 1 , 1 : 1)

        ::

            sage: PP.<t0,t1,t2,t3,t4> = ProductProjectiveSpaces([2,1], ZZ)
            sage: Q = PP([1,1,1,2,1])
            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1,2], ZZ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3+a*b^2, x^2, y^2-z^2, z*y])
            sage: f(Q)
            Traceback (most recent call last):
            ...
            TypeError: (1 : 1 : 1 , 2 : 1) fails to convert into the map's domain
            Product of projective spaces P^1 x P^2 over Integer Ring, but a
            `pushforward` method is not properly implemented
            sage: f([1,1,1,2,1])
            (1 : 2 , 1 : 3 : 2)

        ::

            sage: PP.<x,y,u,v> = ProductProjectiveSpaces(ZZ, [1, 1])
            sage: HP = End(PP)
            sage: g = HP([x^2, y^2, u^2, v^2])
            sage: g([0, 0, 0, 0],check=False)
            (0 : 0 , 0 : 0)
        """
        from sage.schemes.product_projective.point import ProductProjectiveSpaces_point_ring
        if check:
            if not isinstance(P, ProductProjectiveSpaces_point_ring):
                try:
                    P = self.domain()(P)
                except (TypeError, NotImplementedError):
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(P, self.domain()))
            elif self.domain()!= P.codomain():
                raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(P, self.domain()))

        A = self.codomain()
        Q = list(P)
        newP = [f(Q) for f in self.defining_polynomials()]
        return(A.point(newP, check))

    def is_morphism(self):
        r"""
        Returns ``True`` if this mapping is a morphism of products of projective spaces.

        For each component space of the codomain of this mapping we consider the subscheme of
        the domain of this map generated by the corresponding coordinates of the map.
        This map is a morphism if and only if each of these subschemes has no points.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], ZZ)
            sage: H = End(Z)
            sage: f = H([a^2, b^2, x*z-y*z, x^2-y^2, z^2])
            sage: f.is_morphism()
            False

        ::

            sage: P.<x,y,z,u,v,w>=ProductProjectiveSpaces([2, 2], QQ)
            sage: H = End(P)
            sage: f = H([u, v, w, u^2, v^2, w^2])
            sage: f.is_morphism()
            True

        ::

            sage: P.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: Q.<a,b,c,d,e> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H = Hom(P, Q)
            sage: f = H([x^2, y^2, u^3, w^3, u^3])
            sage: f.is_morphism()
            False
        """
        m = 0
        T = self.domain().ambient_space()
        S = self.codomain().ambient_space()

        if T.base_ring().is_field():
            f = self
        else:
            f = self.change_ring(T.base_ring().fraction_field())
            T = T.change_ring(T.base_ring().fraction_field())

        for i in range(S.num_components()):
            t = S[i].dimension_relative() + 1
            X = T.subscheme(list(f)[m : m + t])
            if X.dimension() > -1:
                return False
            m = m + t

        return True

    def nth_iterate(self, P, n, normalize=False):
        r"""
        For a map of this morphism and a point `P` in ``self.domain()``
        this function returns the nth iterate of `P` by this morphism.

        If ``normalize`` is ``True``, then the coordinates are
        automatically normalized.

        .. TODO:: Is there a more efficient way to do this?

        INPUT:

        - ``P`` -- a point in ``self.domain()``.

        - ``n`` -- a positive integer.

        - ``normalize`` - Boolean (optional Default: ``False``).

        OUTPUT:

        - A point in ``self.codomain()``.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y])
            sage: P = Z([1, 1, 1, 1, 1])
            sage: f.nth_iterate(P, 3)
            (1/1872 : 1 , 1 : 1 : 0)
        """
        return(P.nth_iterate(self, n, normalize))

    def orbit(self, P, N, **kwds):
        r"""
        Returns the orbit of `P` by this morphism.

        If `N` is an integer it returns `[P,self(P),\ldots,self^N(P)]`.

        If `n` is a list or tuple `n = [m, k]` it returns `[self^m(P),\ldots,self^k(P)]`.
        Automatically normalize the points if ``normalize == True``. Perform the checks on point initialize if
        ``check==True``.

        INPUT:

        - ``P`` -- a point in ``self.domain()``.

        - ``n`` -- a non-negative integer or list or tuple of two non-negative integers.

        kwds:

        - ``check`` -- boolean (optional - default: ``True``).

        - ``normalize`` -- boolean (optional - default: ``False``).


        OUTPUT:

        - a list of points in ``self.codomain()``.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y])
            sage: P = Z([1, 1, 1, 1, 1])
            sage: f.orbit(P, 3)
            [(1 : 1 , 1 : 1 : 1), (1/2 : 1 , 1 : 0 : 1), (1/12 : 1 , -1 : 1 : 0), (1/1872 : 1 , 1 : 1 : 0)]
        """
        return(P.orbit(self, N, **kwds))

    def nth_iterate_map(self, n):
        r"""
        This function returns the nth iterate of this morphsim as a
        function on ``self.domain()``.

        ALGORITHM:

        Uses a form of successive squaring to reduce computations.


        .. TODO:: This could be improved.

        INPUT:

        - ``n`` -- a positive integer.

        OUTPUT:

        - A map between products of projective spaces.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1 , 2], QQ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3, x^2, y^2, z^2])
            sage: f.nth_iterate_map(3)
            Scheme endomorphism of Product of projective spaces P^1 x P^2 over
            Rational Field
              Defn: Defined by sending (a : b , x : y : z) to
                    (a^27 : b^27 , x^8 : y^8 : z^8).
        """
        if not self.is_endomorphism():
            raise TypeError("domain and codomain of function not equal")

        E = self.domain()
        D = int(n)
        if D < 0:
            raise TypeError("iterate number must be a nonnegative integer")
        N = sum([E.ambient_space()[i].dimension_relative() + 1 for i in range(E.ambient_space().num_components())])
        F = list(self._polys)
        Coord_ring = E.coordinate_ring()
        if isinstance(Coord_ring, QuotientRing_generic):
            PHI = [Coord_ring.gen(i).lift() for i in range(N)]
        else:
            PHI = [Coord_ring.gen(i) for i in range(N)]

        while D:
            if D&1:
                PHI = [PHI[j](*F) for j in range(N)]
            if D > 1: #avoid extra iterate
                F = [F[j](*F) for j in range(N)] #'square'
            D >>= 1
        return End(E)(PHI)
