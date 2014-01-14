r"""
Points on affine varieties

Scheme morphism for points on affine varieties



AUTHORS:

- David Kohel, William Stein

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (2013)
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

from copy                          import copy
from sage.categories.number_fields import NumberFields
_NumberFields = NumberFields()
from sage.rings.integer_ring       import ZZ
from sage.rings.number_field.order import is_NumberFieldOrder
from sage.rings.rational_field     import QQ
from sage.rings.real_mpfr          import RealField
from sage.schemes.generic.morphism import (SchemeMorphism_point, SchemeMorphism, is_SchemeMorphism)
from sage.structure.sequence       import Sequence

############################################################################
# Rational points on schemes, which we view as morphisms determined
# by coordinates.
############################################################################

class SchemeMorphism_point_affine(SchemeMorphism_point):
    """
    A rational point on an affine scheme.

    INPUT:

    - ``X`` -- a subscheme of an ambient affine space over a ring `R`.

    - ``v`` -- a list/tuple/iterable of coordinates in `R`.

    - ``check`` -- boolean (optional, default:``True``). Whether to
      check the input for consistency.

    EXAMPLES::

        sage: A = AffineSpace(2, QQ)
        sage: A(1,2)
        (1, 2)
    """
    def __init__(self, X, v, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_point_affine` for details.

        TESTS::

            sage: from sage.schemes.affine.affine_point import SchemeMorphism_point_affine
            sage: A3.<x,y,z> = AffineSpace(QQ, 3)
            sage: SchemeMorphism_point_affine(A3(QQ), [1,2,3])
            (1, 2, 3)
        """
        SchemeMorphism.__init__(self, X)
        if is_SchemeMorphism(v):
            v = list(v)
        if check:
            # Verify that there are the right number of coords
            d = self.codomain().ambient_space().ngens()
            if len(v) != d:
                raise TypeError("Argument v (=%s) must have %s coordinates."%(v, d))
            if not isinstance(v,(list,tuple)):
                raise TypeError("Argument v (= %s) must be a scheme point, list, or tuple."%str(v))
            # Make sure the coordinates all lie in the appropriate ring
            v = Sequence(v, X.value_ring())
            # Verify that the point satisfies the equations of X.
            X.extended_codomain()._check_satisfies_equations(v)
        self._coords = tuple(v)

    def nth_iterate(self,f,n):
        r"""
        Returns the point `f^n(self)`

        INPUT:

        - ``f`` -- a :class:`SchemeMorphism_polynomial` with ``self`` if ``f.domain()``
        - ``n`` -- a positive integer.

        OUTPUT:

        - a point in ``f.codomain()``

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(QQ,2)
            sage: H=Hom(A,A)
            sage: f=H([(x-2*y^2)/x,3*x*y])
            sage: A(9,3).nth_iterate(f,3)
            (-104975/13123, -9566667)

        ::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: X=A.subscheme([x-y^2])
            sage: H=Hom(X,X)
            sage: f=H([9*y^2,3*y])
            sage: X(9,3).nth_iterate(f,4)
            (59049, 243)
        """
        if self.codomain()!=f.domain():
            raise TypeError("Point is not defined over domain of function")
        if f.domain() != f.codomain():
            raise TypeError("Domain and Codomain of function not equal")
        if n==0:
            return(self)
        else:
            Q=f(self)
            for i in range(2,n+1):
                Q=f(Q)
            return(Q)

    def orbit(self,f,N):
        r"""
        Returns the orbit of self by `f`. If `n` is an integer it returns `[self,f(self),\ldots,f^{n}(self)]`.

        If `n` is a list or tuple `n=[m,k]` it returns `[f^{m}(self),\ldots,f^{k}(self)]`.

        INPUT:

        - ``f`` -- a :class:`SchemeMorphism_polynomial` with ``self`` in ``f.domain()``
        - ``n`` -- a non-negative integer or list or tuple of two non-negative integers

        OUTPUT:

        - a list of points in ``f.codomain()``

        EXAMPLES::

            sage: A.<x,y>=AffineSpace(QQ,2)
            sage: H=Hom(A,A)
            sage: f=H([(x-2*y^2)/x,3*x*y])
            sage: A(9,3).orbit(f,3)
            [(9, 3), (-1, 81), (13123, -243), (-104975/13123, -9566667)]

        ::

            sage: A.<x>=AffineSpace(QQ,1)
            sage: H=Hom(A,A)
            sage: f=H([(x-2)/x])
            sage: A(1/2).orbit(f,[1,3])
            [(-3), (5/3), (-1/5)]

        ::

            sage: A.<x,y>=AffineSpace(ZZ,2)
            sage: X=A.subscheme([x-y^2])
            sage: H=Hom(X,X)
            sage: f=H([9*y^2,3*y])
            sage: X(9,3).orbit(f,(0,4))
            [(9, 3), (81, 9), (729, 27), (6561, 81), (59049, 243)]
        """
        Q=copy(self)
        if type(N)==list or type(N)==tuple:
            Bounds=list(N)
        else:
            Bounds=[0,N]
        for i in range(1,Bounds[0]+1):
            Q=f(Q)
        Orb=[Q]
        for i in range(Bounds[0]+1,Bounds[1]+1):
            Q=f(Q)
            Orb.append(Q)
        return(Orb)

    def global_height(self, prec=None):
        r"""
        Returns the logarithmic height of the point.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y>=AffineSpace(QQ,2)
            sage: Q=P(41,1/12)
            sage: Q.global_height()
            3.71357206670431

        ::

            sage: P=AffineSpace(ZZ,4,'x')
            sage: Q=P(3,17,-51,5)
            sage: Q.global_height()
            3.93182563272433

        ::

            sage: R.<x>=PolynomialRing(QQ)
            sage: k.<w>=NumberField(x^2+5)
            sage: A=AffineSpace(k,2,'z')
            sage: A([3,5*w+1]).global_height(prec=100)
            2.4181409534757389986565376694

        .. TODO::

            p-adic heights

            add heights to integer.pyx and remove special case
        """
        if self.domain().base_ring() == ZZ:
            if prec is None:
                R = RealField()
            else:
                R = RealField(prec)
            H=max([self[i].abs() for i in range(self.codomain().ambient_space().dimension_relative())])
            return(R(max(H,1)).log())
        if self.domain().base_ring() in _NumberFields or is_NumberFieldOrder(self.domain().base_ring()):
            return(max([self[i].global_height(prec) for i in range(self.codomain().ambient_space().dimension_relative())]))
        else:
            raise NotImplementedError("Must be over a Numberfield or a Numberfield Order")

class SchemeMorphism_point_affine_field(SchemeMorphism_point_affine):
    pass

class SchemeMorphism_point_affine_finite_field(SchemeMorphism_point_affine_field):

    def __hash__(self):
        r"""
        Returns the integer hash of ``self``

        OUTPUT:

        - integer

        EXAMPLES::

            sage: P.<x,y,z>=AffineSpace(GF(5),3)
            sage: hash(P(2,1,2))
            57

        ::

            sage: P.<x,y,z>=AffineSpace(GF(7),3)
            sage: X=P.subscheme(x^2-y^2)
            sage: hash(X(1,1,2))
            106

        ::

            sage: P.<x,y>=AffineSpace(GF(13),2)
            sage: hash(P(3,4))
            55

        ::

            sage: P.<x,y>=AffineSpace(GF(13^3,'t'),2)
            sage: hash(P(3,4))
            8791
        """
        p=self.codomain().base_ring().order()
        N=self.codomain().ambient_space().dimension_relative()
        return sum(hash(self[i])*p**i for i in range(N))

    def orbit_structure(self,f):
        r"""
        Every points is preperiodic over a finite field. This funtion returns the pair `[m,n]` where `m` is the
        preperiod and `n` the period of the point ``self`` by ``f``.

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        OUTPUT:

        - a list `[m,n]` of integers

        EXAMPLES::

            sage: P.<x,y,z>=AffineSpace(GF(5),3)
            sage: H=Hom(P,P)
            sage: f=H([x^2+y^2,y^2,z^2+y*z])
            sage: P(1,1,1).orbit_structure(f)
            [0, 6]

        ::

            sage: P.<x,y,z>=AffineSpace(GF(7),3)
            sage: X=P.subscheme(x^2-y^2)
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,z^2])
            sage: X(1,1,2).orbit_structure(f)
            [0, 2]

        ::

            sage: P.<x,y>=AffineSpace(GF(13),2)
            sage: H=Hom(P,P)
            sage: f=H([x^2-y^2,y^2])
            sage: P(3,4).orbit_structure(f)
            [2, 6]
        """
        Orbit=[]
        index=1
        P=copy(self)
        F=copy(f)
        while not P in Orbit:
            Orbit.append(P)
            P=F(P)
            index+=1
        I=Orbit.index(P)
        return([I,index-I-1])

