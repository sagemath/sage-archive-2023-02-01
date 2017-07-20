r"""
Dynamical systems for products of projective spaces

This class builds on the prouct projective space class.

EXAMPLES::

    sage: P1xP1.<x,y,u,v> = ProductProjectiveSpaces(QQ, [1, 1])
    sage: DynamicalSystem_projective([x^2*u, y^2*v, x*v^2, y*u^2], domain=P1xP1)
    Dynamical System of Product of projective spaces P^1 x P^1 over Rational Field
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
from copy import copy
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem_generic
from sage.rings.all import ZZ
from sage.rings.quotient_ring import QuotientRing_generic
from sage.schemes.product_projective.morphism import ProductProjectiveSpaces_morphism_ring


class DynamicalSystem_product_projective_ring(DynamicalSystem_generic,\
                                             ProductProjectiveSpaces_morphism_ring):
    r"""
    The class of dynamical systems on products of projective spaces.

    The components are projective space dynamical systems.

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: DynamicalSystem_projective([x^2, y^2, z^2, w^2, u^2], domain=T)
        Dynamical System of Product of projective spaces P^2 x P^1 over Rational Field
              Defn: Defined by sending (x : y : z , w : u) to
                    (x^2 : y^2 : z^2 , w^2 : u^2).
    """

    def __init__(self, polys, domain, check=True):
        r"""
        The Python constructor.

        INPUT:

        - ``polys`` -- anything that defines a point in the class.

        - ``domain`` -- product of projective spaces scheme or subscheme

        - ``check`` -- Boolean. Whether or not to perform input checks.
          (Default:`` True``)

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u^2], domain=T)
            Dynamical System of Product of projective spaces P^2 x P^1 over Rational Field
              Defn: Defined by sending (x : y : z , w : u) to
                    (x^2*u : y^2*w : z^2*u , w^2 : u^2).

        ::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u*z], domain=T)
            Traceback (most recent call last):
            ...
            TypeError: polys (=[x^2*u, y^2*w, z^2*u, w^2, z*u]) must be
            multi-homogeneous of the same degrees (by component)
        """
        if check:
            #check multi-homogeneous
            #if self is a subscheme, we may need the lift of the polynomials
            try:
                polys[0].exponents()
            except AttributeError:
                polys = [f.lift() for f in polys]

            splitpolys = domain._factors(polys)
            for m in range(len(splitpolys)):
                d = domain._degree(splitpolys[m][0])
                if not all(d == domain._degree(f) for f in splitpolys[m]):
                    raise  TypeError("polys (=%s) must be multi-homogeneous of the same degrees (by component)"%polys)

        DynamicalSystem_generic.__init__(self, polys, domain=domain, check=check)

    def __call__(self, P, check=True):
        r"""
        Make dynamical systemsof products of projective spaces callable.

        INPUT:

        - ``P`` -- a point in the domain.

        - ``check`` -- Boolean - whether or not to perform the input checks
          on the image point (Default: ``True``).

        OUTPUT: The image point in the codomain.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: F = DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u^2], domain=T)
            sage: F(T([2, 1, 3, 0, 1]))
            (4/9 : 0 : 1 , 0 : 1)
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

        A = self.domain()
        Q = list(P)
        newP = [f(Q) for f in self.defining_polynomials()]
        return(A.point(newP, check))

    def nth_iterate(self, P, n, normalize=False):
        r"""
        Return the ``n``-th iterate of ``P`` by this dynamical system.

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
            sage: f = DynamicalSystem_projective([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y], domain=Z)
            sage: P = Z([1, 1, 1, 1, 1])
            sage: f.nth_iterate(P, 3)
            (1/1872 : 1 , 1 : 1 : 0)

        ::

            sage: Z.<a,b,x,y> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: f = DynamicalSystem_projective([a*b, b^2, x^3 - y^3, y^2*x], domain=Z)
            sage: P = Z([2, 6, 2, 4])
            sage: f.nth_iterate(P, 2, normalize = True)
            (1 : 3 , 407 : 112)
        """
        if P.codomain() != self.domain():
            raise TypeError("point is not defined over domain of function")
        n = ZZ(n)
        if n < 0:
            raise TypeError("must be a forward orbit")
        if n == 0:
            return(self)
        else:
            Q = self(P)
            if normalize == True:
                Q.normalize_coordinates()
            for i in range(2,n+1):
                Q = self(Q)
                if normalize == True:
                    Q.normalize_coordinates()
            return(Q)

    def orbit(self, P, N, **kwds):
        r"""
        Return the orbit of `P` by this dynamcial system.

        If `N` is an integer it returns `[P,self(P),\ldots,self^N(P)]`.

        If `N` is a list or tuple `N = [m, k]` it returns `[self^m(P),\ldots,self^k(P)]`.
        Automatically normalize the points if ``normalize == True``. Perform the checks on point initialize if
        ``check==True``.

        INPUT:

        - ``P`` -- a point in ``self.domain()``.

        - ``N`` -- a non-negative integer or list or tuple of two non-negative integers.

        kwds:

        - ``check`` -- boolean (optional - default: ``True``).

        - ``normalize`` -- boolean (optional - default: ``False``).


        OUTPUT:

        - a list of points in ``self.codomain()``.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: f = DynamicalSystem_projective([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y], domain=Z)
            sage: P = Z([1, 1, 1, 1, 1])
            sage: f.orbit(P, 3)
            [(1 : 1 , 1 : 1 : 1), (1/2 : 1 , 1 : 0 : 1), (1/12 : 1 , -1 : 1 : 0), (1/1872 : 1 , 1 : 1 : 0)]

        ::

            sage: Z.<a,b,x,y> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: f = DynamicalSystem_projective([a*b, b^2, x^3 - y^3, y^2*x], domain=Z)
            sage: P = Z([2, 6, 2, 4])
            sage: f.orbit(P, 3, normalize=True)
            [(1 : 3 , 1 : 2), (1 : 3 , -7 : 4), (1 : 3 , 407 : 112), (1 : 3 , 66014215 : 5105408)]
        """
        if P.codomain() != self.domain():
            raise TypeError("point is not defined over domain of function")
        if not isinstance(N, (list,tuple)):
            N = [0, N]
        try:
            N[0] = ZZ(N[0])
            N[1] = ZZ(N[1])
        except TypeError:
            raise TypeError("orbit bounds must be integers")
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be non-negative")
        if N[0] > N[1]:
            return([])

        Q = copy(P)
        check = kwds.pop("check", True)
        normalize = kwds.pop("normalize", False)

        if normalize == True:
            Q.normalize_coordinates()
        for i in range(1, N[0]+1):
            Q = self(Q, check)
            if normalize == True:
                Q.normalize_coordinates()
        orb = [Q]
        for i in range(N[0]+1, N[1]+1):
            Q = self(Q, check)
            if normalize == True:
                Q.normalize_coordinates()
            orb.append(Q)
        return(orb)

    def nth_iterate_map(self, n):
        r"""
        Return the nth iterate of this dynamical system.

        ALGORITHM:

        Uses a form of successive squaring to reduce computations.


        .. TODO:: This could be improved.

        INPUT:

        - ``n`` -- a positive integer.

        OUTPUT:

        - A dynamical system of products of projective spaces.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1 , 2], QQ)
            sage: f = DynamicalSystem_projective([a^3, b^3, x^2, y^2, z^2], domain=Z)
            sage: f.nth_iterate_map(3)
            Dynamical System of Product of projective spaces P^1 x P^2 over
            Rational Field
              Defn: Defined by sending (a : b , x : y : z) to
                    (a^27 : b^27 , x^8 : y^8 : z^8).
        """
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
        return DynamicalSystem_product_projective_ring(PHI, domain=self.domain())
