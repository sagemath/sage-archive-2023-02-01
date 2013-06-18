r"""
Points on projective varieties

Scheme morphism for points on projective varieties



AUTHORS:

- David Kohel, William Stein

- William Stein (2006-02-11): fixed bug where P(0,0,0) was allowed as
  a projective point.

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz: (June 2012) added support for projective ring
            (March 2013) iteration functionality and new directory structure
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


from sage.rings.infinity       import infinity
from sage.rings.arith          import gcd, lcm
from sage.rings.integer_ring   import ZZ
from sage.rings.quotient_ring  import QuotientRing_generic
from copy                      import copy
from sage.schemes.generic.morphism import (SchemeMorphism,
                                           is_SchemeMorphism,
                                           SchemeMorphism_point)
from sage.structure.element    import AdditiveGroupElement
from sage.structure.sequence   import Sequence



#*******************************************************************
# Projective varieties
#*******************************************************************
class SchemeMorphism_point_projective_ring(SchemeMorphism_point):
    """
    A rational point of projective space over a ring.

    INPUT:

    -  ``X`` -- a homset of a subscheme of an ambient projective space over a field `K`

    - ``v`` -- a list or tuple of coordinates in `K`

    - ``check`` -- boolean (optional, default:``True``). Whether to check the input for consistency.

    EXAMPLES::

        sage: P = ProjectiveSpace(2, ZZ)
        sage: P(2,3,4)
        (2 : 3 : 4)

    """

    def __init__(self, X, v, check=True):
        """
        The Python constructor.

        EXAMPLES::

            sage: P = ProjectiveSpace(2, QQ)
            sage: P(2, 3/5, 4)
            (1/2 : 3/20 : 1)

        ::

            sage: P = ProjectiveSpace(1, ZZ)
            sage: P([0, 1])
            (0 : 1)

        ::

            sage: P = ProjectiveSpace(1, ZZ)
            sage: P([0, 0, 1])
            Traceback (most recent call last):
            ...
            TypeError: v (=[0, 0, 1]) must have 2 components

        ::

            sage: P = ProjectiveSpace(3, QQ)
            sage: P(0,0,0,0)
            Traceback (most recent call last):
            ...
            ValueError: [0, 0, 0, 0] does not define a valid point since all entries are 0

        ::

        It is possible to avoid the possibly time consuming checks, but be careful!!

            sage: P = ProjectiveSpace(3, QQ)
            sage: P.point([0,0,0,0],check=False)
            (0 : 0 : 0 : 0)

        ::

            sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
            sage: X=P.subscheme([x^2-y*z])
            sage: X([2,2,2])
            (2 : 2 : 2)
        """
        SchemeMorphism.__init__(self, X)
        if check:
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            d = X.codomain().ambient_space().ngens()
            if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            elif v is infinity:
                v = [0] * (d)
                v[1] = 1
            if not isinstance(v,(list,tuple)):
                raise TypeError("Argument v (= %s) must be a scheme point, list, or tuple."%str(v))
            if len(v) != d and len(v) != d-1:
                raise TypeError("v (=%s) must have %s components"%(v, d))

            R = X.value_ring()
            v = Sequence(v, R)
            if len(v) == d-1:     # very common special case
                v.append(1)

            n = len(v)
            all_zero = True
            for i in range(n):
                last = n-1-i
                if v[last]:
                    all_zero = False
                    break
            if all_zero:
                raise ValueError("%s does not define a valid point since all entries are 0"%repr(v))

            X.extended_codomain()._check_satisfies_equations(v)

            if isinstance(X.codomain().base_ring(), QuotientRing_generic):
                lift_coords = [P.lift() for P in v]
            else:
                lift_coords = v
            v = Sequence(lift_coords)

        self._coords = v

    def __eq__(self,right):
        """
        Tests the proejctive equality of two points.

        INPUT:

        - ``right`` - a point on projective space

        OUTPUT:

        - Boolean - True if ``self`` and ``right`` define the same point. False otherwise.

        Examples::

            sage: PS=ProjectiveSpace(ZZ,1,'x')
            sage: P=PS([1,2])
            sage: Q=PS([2,4])
            sage: P==Q
            True

        ::

            sage: PS=ProjectiveSpace(ZZ,1,'x')
            sage: P=PS([1,2])
            sage: Q=PS([1,0])
            sage: P==Q
            False

        ::

            sage: PS=ProjectiveSpace(Zp(5),1,'x')
            sage: P=PS([0,1])
            sage: P==0
            True

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: PS=ProjectiveSpace(R,1,'x')
            sage: P=PS([t,1+t^2])
            sage: Q=PS([t^2, t+t^3])
            sage: P==Q
            True

        ::

            sage: PS=ProjectiveSpace(ZZ,2,'x')
            sage: P=PS([0,1,2])
            sage: P==0
            False

        ::

            sage: PS=ProjectiveSpace(ZZ,1,'x')
            sage: P=PS([2,1])
            sage: PS2=ProjectiveSpace(Zp(7),1,'x')
            sage: Q=PS2([2,1])
            sage: P==Q
            False

        ::

            sage: PS=ProjectiveSpace(ZZ.quo(6),2,'x')
            sage: P=PS([2,4,1])
            sage: Q=PS([0,1,3])
            sage: P==Q
            False

        ::

            sage: R.<z>=PolynomialRing(QQ)
            sage: K.<t>=NumberField(z^2+5)
            sage: OK=K.ring_of_integers()
            sage: t=OK.gen(1)
            sage: PS.<x,y>=ProjectiveSpace(OK,1)
            sage: P=PS(2,1+t)
            sage: Q=PS(1-t,3)
            sage: P==Q
            True

        """
        if not isinstance(right, SchemeMorphism_point):
            try:
                right = self.codomain()(right)
            except TypeError:
                return False
        if self.codomain()!=right.codomain():
            return False
        n=len(self._coords)
        for i in range(0,n):
            for j in range(i+1,n):
                if self._coords[i]*right._coords[j] != self._coords[j]*right._coords[i]:
                    return False
        return True

    def __ne__(self,right):
        """
        Tests the proejctive equality of two points.

        INPUT:

        - ``right`` - a point on projective space

        OUTPUT:

        - Boolean - True if ``self`` and ``right`` define different points. False otherwise.

        Examples::

            sage: PS=ProjectiveSpace(ZZ,1,'x')
            sage: P=PS([1,2])
            sage: Q=PS([2,4])
            sage: P!=Q
            False

        ::

            sage: PS=ProjectiveSpace(ZZ,1,'x')
            sage: P=PS([1,2])
            sage: Q=PS([1,0])
            sage: P!=Q
            True

        ::

            sage: PS=ProjectiveSpace(Zp(5),1,'x')
            sage: P=PS([0,1])
            sage: P!=0
            False

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: PS=ProjectiveSpace(R,1,'x')
            sage: P=PS([t,1+t^2])
            sage: Q=PS([t^2, t+t^3])
            sage: P!=Q
            False

        ::

            sage: PS=ProjectiveSpace(ZZ,2,'x')
            sage: P=PS([0,1,2])
            sage: P!=0
            True

        ::

            sage: PS=ProjectiveSpace(ZZ,1,'x')
            sage: P=PS([2,1])
            sage: PS2=ProjectiveSpace(Zp(7),1,'x')
            sage: Q=PS2([2,1])
            sage: P!=Q
            True

        ::

            sage: PS=ProjectiveSpace(ZZ.quo(6),2,'x')
            sage: P=PS([2,4,1])
            sage: Q=PS([0,1,3])
            sage: P!=Q
            True
        """
        if not isinstance(right, SchemeMorphism_point):
            try:
                right = self.codomain()(right)
            except TypeError:
                return True
        if self.codomain()!=right.codomain():
            return True
        n=len(self._coords)
        for i in range(0,n):
            for j in range(i+1,n):
                if self._coords[i]*right._coords[j] != self._coords[j]*right._coords[i]:
                    return True
        return False

    def scale_by(self,t):
        """
        Scale the coordinates of the point ``self`` by `t`. A ``TypeError`` occurs if
        the point is not in the base_ring of the codomain after scaling.

        INPUT:

        - ``t`` -- a ring element

        OUTPUT:

        - None.

        EXAMPLES::

            sage: R.<t>=PolynomialRing(QQ)
            sage: P=ProjectiveSpace(R,2,'x')
            sage: p=P([3/5*t^3,6*t, t])
            sage: p.scale_by(1/t); p
            (3/5*t^2 : 6 : 1)

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: S=R.quo(R.ideal(t^3))
            sage: P.<x,y,z>=ProjectiveSpace(S,2)
            sage: Q=P(t,1,1)
            sage: Q.scale_by(t);Q
            (tbar^2 : tbar : tbar)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: Q=P(2,2,2)
            sage: Q.scale_by(1/2);Q
            (1 : 1 : 1)
        """
        if t==0:  #what if R(t) == 0 ?
            raise ValueError("Cannot scale by 0")
        R=self.codomain().base_ring()
        if isinstance(R, QuotientRing_generic):
            phi=R.coerce_map_from(self.codomain().ambient_space().base_ring())
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                self._coords[i]=phi(self._coords[i]*t).lift()
        else:
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                self._coords[i]=R(self._coords[i]*t)

    def normalize_coordinates(self):
        """
        Removes common factors from the coordinates of ``self`` (including `-1`).

        OUTPUT:

        - None.

        EXAMPLES::

            sage: P=ProjectiveSpace(ZZ,2,'x')
            sage: p=P([-5,-15,-20])
            sage: p.normalize_coordinates(); p
            (1 : 3 : 4)

        ::

            sage: P=ProjectiveSpace(Zp(7),2,'x')
            sage: p=P([-5,-15,-2])
            sage: p.normalize_coordinates(); p
            (5 + O(7^20) : 1 + 2*7 + O(7^20) : 2 + O(7^20))

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: P=ProjectiveSpace(R,2,'x')
            sage: p=P([3/5*t^3,6*t, t])
            sage: p.normalize_coordinates(); p
            (3/5*t^2 : 6 : 1)

        ::

            sage: P.<x,y>=ProjectiveSpace(Zmod(20),1)
            sage: Q=P(4,8)
            sage: Q.normalize_coordinates()
            sage: Q
            (1 : 2)

        ::

            sage: R.<t>=PolynomialRing(QQ,1)
            sage: S=R.quotient_ring(R.ideal(t^3))
            sage: P.<x,y>=ProjectiveSpace(S,1)
            sage: Q=P(t,t^2)
            sage: Q.normalize_coordinates()
            sage: Q
            (1 : t)
        """
        R=self.codomain().base_ring()
        GCD = R(gcd(self[0],self[1]))
        index=2
        if self[0]>0 or self[1] >0:
            neg=0
        else:
            neg=1
        while GCD!=1 and index < len(self._coords):
            if self[index]>0:
                neg=0
            GCD=R(gcd(GCD,self[index]))
            index+=1
        if isinstance(R,(QuotientRing_generic)):
            R=R.cover_ring()
            GCD=GCD.lift()
        if GCD != 1:
            if neg==1:
                self.scale_by(R(-1)/GCD)
            else:
                self.scale_by(R(1)/GCD)
        else:
            if neg==1:
                self.scale_by(R(-1))

    def dehomogenize(self,n):
        r"""
        Dehomogenizes at the nth coordinate

        INPUT:

        - ``n`` -- non-negative integer

        OUTPUT:

        - :class:`SchemeMorphism_point_affine`

        EXAMPLES::

            sage: P.<x,y,z>=ProjectiveSpace(QQ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: Q=X(23,23,46)
            sage: Q.dehomogenize(2)
            (1/2, 1/2)

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: S=R.quo(R.ideal(t^3))
            sage: P.<x,y,z>=ProjectiveSpace(S,2)
            sage: Q=P(t,1,1)
            sage: Q.dehomogenize(1)
            (tbar, 1)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(GF(5),2)
            sage: Q=P(1,3,1)
            sage: Q.dehomogenize(0)
            (3, 1)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(GF(5),2)
            sage: Q=P(1,3,0)
            sage: Q.dehomogenize(2)
            Traceback (most recent call last):
            ...
            ValueError: Can't dehomogenize at 0 coordinate.
        """
        if self[n]==0:
            raise ValueError("Can't dehomogenize at 0 coordinate.")
        PS=self.codomain()
        A=PS.affine_patch(n)
        Q=[]
        for i in range(0,PS.ambient_space().dimension_relative()+1):
            if i !=n:
                Q.append(self[i]/self[n])
        return(A.point(Q))

    def nth_iterate(self,f,n,normalize=False):
        r"""
        For a map ``self`` and a point `P` in ``self.domain()``
        this function returns the nth iterate of `P` by ``self``. If ``normalize==True``,
        then the coordinates are automatically normalized.

        INPUT:

        - ``f`` -- a SchmemMorphism_polynomial with ``self`` in ``f.domain()``

        - ``n`` -- a positive integer.

        - ``normalize`` -- Boolean (optional Default: ``False``)

        OUTPUT:

        - A point in ``self.codomain()``

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(ZZ,1)
            sage: H=Hom(P,P)
            sage: f=H([x^2+y^2,2*y^2])
            sage: P(1,1).nth_iterate(f,4)
            (32768 : 32768)

        ::

            sage: P.<x,y>=ProjectiveSpace(ZZ,1)
            sage: H=Hom(P,P)
            sage: f=H([x^2+y^2,2*y^2])
            sage: P(1,1).nth_iterate(f,4,1)
            (1 : 1)

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: P.<x,y,z>=ProjectiveSpace(R,2)
            sage: H=Hom(P,P)
            sage: f=H([x^2+t*y^2,(2-t)*y^2,z^2])
            sage: P(2+t,7,t).nth_iterate(f,2)
            (t^4 + 2507*t^3 - 6787*t^2 + 10028*t + 16 : -2401*t^3 + 14406*t^2 -
            28812*t + 19208 : t^4)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: X=P.subscheme(x^2-y^2)
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,z^2])
            sage: X(2,2,3).nth_iterate(f,3)
            (256 : 256 : 6561)

        .. TODO:: Is there a more efficient way to do this?
        """
        if self.codomain()!=f.domain():
            raise TypeError("Point is not defined over domain of function")
        if f.domain() != f.codomain():
            raise TypeError("Domain and Codomain of function not equal")
        try:
            n=ZZ(n)
        except TypeError:
            raise TypeError("Iterate number must be an integer")
        if n <0:
            raise TypeError("Must be a forward orbit")
        if n==0:
            return(self)
        else:
            Q=f(self)
            if normalize==True:
                Q.normalize_coordinates()
            for i in range(2,n+1):
                Q=f(Q)
                if normalize==True:
                    Q.normalize_coordinates()
            return(Q)

    def orbit(self,f,N,**kwds):
        r"""
        Returns the orbit of `P` by ``self``. If `n` is an integer it returns `[P,self(P),\ldots,self^n(P)]`.
        If `n` is a list or tuple `n=[m,k]` it returns `[self^m(P),\ldots,self^k(P)`].
        Automatically normalize the points if ``normalize=True``. Perform the checks on point initialization if
        ``check=True``

        INPUT:

        - ``f`` -- a :class:`SchemeMorphism_polynomial` with ``self`` in ``f.domain()``

        - ``N`` -- a non-negative integer or list or tuple of two non-negative integers

        kwds:

        - ``check`` -- boolean (optional - default: ``True``)

        - ``normalize`` -- boolean (optional - default: ``False``)


        OUTPUT:

        - a list of points in ``self.codomain()``

        EXAMPLES::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: H=Hom(P,P)
            sage: f=H([x^2+y^2,y^2-z^2,2*z^2])
            sage: P(1,2,1).orbit(f,3)
            [(1 : 2 : 1), (5 : 3 : 2), (34 : 5 : 8), (1181 : -39 : 128)]

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: H=Hom(P,P)
            sage: f=H([x^2+y^2,y^2-z^2,2*z^2])
            sage: P(1,2,1).orbit(f,[2,4])
            [(34 : 5 : 8), (1181 : -39 : 128), (1396282 : -14863 : 32768)]

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: X=P.subscheme(x^2-y^2)
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,x*z])
            sage: X(2,2,3).orbit(f,3,normalize=True)
            [(2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3)]

        ::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(P,P)
            sage: f=H([x^2+y^2,y^2])
            sage: P.point([1,2],False).orbit(f,4,check=False)
            [(1 : 2), (5 : 4), (41 : 16), (1937 : 256), (3817505 : 65536)]
        """
        if (isinstance(N,(list,tuple))==False):
            N=[0,N]
        try:
            N[0]=ZZ(N[0])
            N[1]=ZZ(N[1])
        except TypeError:
            raise TypeError("Orbit bounds must be integers")
        if N[0]<0 or N[1] <0:
            raise TypeError("Orbit bounds must be non-negative")
        if N[0] > N[1]:
            return([])

        Q=copy(self)
        check = kwds.pop("check",True)
        normalize = kwds.pop("normalize",False)

        if normalize==True:
            Q.normalize_coordinates()
        for i in range(1,N[0]+1):
            Q=f(Q,check)
            if normalize==True:
                Q.normalize_coordinates()
        Orb=[Q]
        for i in range(N[0]+1,N[1]+1):
            Q=f(Q,check)
            if normalize==True:
                Q.normalize_coordinates()
            Orb.append(Q)
        return(Orb)


class SchemeMorphism_point_projective_field(SchemeMorphism_point_projective_ring):
    """
    A rational point of projective space over a field.

    INPUT:

    -  ``X`` -- a homset of a subscheme of an ambient projective space
       over a field `K`

    - ``v`` -- a list or tuple of coordinates in `K`

    - ``check`` -- boolean (optional, default:``True``). Whether to
      check the input for consistency.

    EXAMPLES::

        sage: P = ProjectiveSpace(3, RR)
        sage: P(2,3,4,5)
        (0.400000000000000 : 0.600000000000000 : 0.800000000000000 : 1.00000000000000)
    """

    def __init__(self, X, v, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_point_projective_ring` for details.

        This function still normalized points so that the rightmost non-zero coordinate is 1. The is to maintain current functionality with current
        implementations of curves in projectives space (plane, connic, elliptic, etc). The class:`SchemeMorphism_point_projective_ring` is for general use.

        EXAMPLES::

            sage: P = ProjectiveSpace(2, QQ)
            sage: P(2, 3/5, 4)
            (1/2 : 3/20 : 1)

        ::

            sage: P = ProjectiveSpace(3, QQ)
            sage: P(0,0,0,0)
            Traceback (most recent call last):
            ...
            ValueError: [0, 0, 0, 0] does not define a valid point since all entries are 0

        ::

            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: X=P.subscheme([x^2-y*z])
            sage: X([2,2,2])
            (1 : 1 : 1)
        """
        SchemeMorphism.__init__(self, X)
        if check:
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            d = X.codomain().ambient_space().ngens()
            if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            elif v is infinity:
                v = [0] * (d)
                v[1] = 1
            if not isinstance(v,(list,tuple)):
                raise TypeError("Argument v (= %s) must be a scheme point, list, or tuple."%str(v))
            if len(v) != d and len(v) != d-1:
                raise TypeError("v (=%s) must have %s components"%(v, d))

            R = X.value_ring()
            v = Sequence(v, R)
            if len(v) == d-1:     # very common special case
                v.append(1)

            n = len(v)
            all_zero = True
            for i in range(n):
                last = n-1-i
                if v[last]:
                    all_zero = False
                    c = v[last]
                    if c == R.one():
                        break
                    for j in range(last):
                        v[j] /= c
                    v[last] = R.one()
                    break
            if all_zero:
                raise ValueError("%s does not define a valid point since all entries are 0"%repr(v))

            X.extended_codomain()._check_satisfies_equations(v)

        self._coords = v

    def normalize_coordinates(self):
        r"""
        Normalizes ``self`` so that the last non-zero coordinate is `1`.

        OUTPUT: None.

        EXAMPLES::

            sage: P.<x,y,z>=ProjectiveSpace(GF(5),2)
            sage: Q=P.point([1,3,0],false);Q
            (1 : 3 : 0)
            sage: Q.normalize_coordinates();Q
            (2 : 1 : 0)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(QQ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: Q=X.point([23,23,46], false);Q
            (23 : 23 : 46)
            sage: Q.normalize_coordinates();Q
            (1/2 : 1/2 : 1)
        """
        index=self.codomain().ambient_space().dimension_relative()
        while self[index]==0:
            index-=1
        self.scale_by(1/self[index])


    def clear_denominators(self):
        r"""
        scales by the least common multiple of the denominators.

        OUTPUT: None.

        EXAMPLES::

            sage: R.<t>=PolynomialRing(QQ)
            sage: P.<x,y,z>=ProjectiveSpace(FractionField(R),2)
            sage: Q=P([t,3/t^2,1])
            sage: Q.clear_denominators(); Q
            (t^3 : 3 : t^2)

        ::

            sage: R.<x>=PolynomialRing(QQ)
            sage: K.<w>=NumberField(x^2-3)
            sage: P.<x,y,z>=ProjectiveSpace(K,2)
            sage: Q=P([1/w,3,0])
            sage: Q.clear_denominators(); Q
            (w : 9 : 0)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(QQ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: Q=X([1/2,1/2,1]);
            sage: Q.clear_denominators(); Q
            (1 : 1 : 2)
        """
        self.scale_by(lcm([self[i].denominator() for i in range(self.codomain().ambient_space().dimension_relative())]))

class SchemeMorphism_point_projective_finite_field(SchemeMorphism_point_projective_field):

    def __hash__(self):
        r"""
        Returns the integer hash of ``self``


        OUTPUT:

        - integer

        EXAMPLES::

            sage: P.<x,y,z>=ProjectiveSpace(GF(5),2)
            sage: hash(P(2,1,2))
            41

        ::

            sage: P.<x,y,z>=ProjectiveSpace(GF(7),2)
            sage: X=P.subscheme(x^2-y^2)
            sage: hash(X(1,1,2))
            81

        ::

            sage: P.<x,y>=ProjectiveSpace(GF(13),1)
            sage: hash(P(3,4))
            17

        ::

            sage: P.<x,y>=ProjectiveSpace(GF(13^3,'t'),1)
            sage: hash(P(3,4))
            2201
        """
        p=self.codomain().base_ring().order()
        N=self.codomain().ambient_space().dimension_relative()
        return sum(hash(self[i])*p**i for i in range(N+1))

    def orbit_structure(self,f):
        r"""
        Every points is preperiodic over a finite field. This funtion returns the pair `[m,n]` where `m` is the
        preperiod and `n` the period of the point `P` by ``self``.

        INPUT:

        - ``f`` -- a :class:`ScemeMorphism_polynomial` with ``self`` in ``f.domain()``

        OUTPUT:

        - a list `[m,n]` of integers

        EXAMPLES::

            sage: P.<x,y,z>=ProjectiveSpace(GF(5),2)
            sage: H=Hom(P,P)
            sage: f=H([x^2+y^2,y^2,z^2 + y*z])
            sage: P(1,0,1).orbit_structure(f)
            [0, 1]

        ::

            sage: P.<x,y,z>=ProjectiveSpace(GF(17),2)
            sage: X=P.subscheme(x^2-y^2)
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,z^2])
            sage: X(1,1,2).orbit_structure(f)
            [3, 1]

        ::

            sage: R.<t> = GF(13^3)
            sage: P.<x,y>=ProjectiveSpace(R,1)
            sage: H=Hom(P,P)
            sage: f=H([x^2-y^2,y^2])
            sage: P(t,4).orbit_structure(f)
            [11, 6]
        """
        Orbit=[]
        index=1
        P=copy(self)
        P.normalize_coordinates()
        F=copy(f)
        F.normalize_coordinates()
        while not P in Orbit:
            Orbit.append(P)
            P=F(P)
            P.normalize_coordinates()
            index+=1
        I=Orbit.index(P)
        return([I,index-I-1])

#*******************************************************************
# Abelian varieties
#*******************************************************************
class SchemeMorphism_point_abelian_variety_field(AdditiveGroupElement, SchemeMorphism_point_projective_field):
    """
    A rational point of an abelian variety over a field.

    EXAMPLES::

        sage: E = EllipticCurve([0,0,1,-1,0])
        sage: origin = E(0)
        sage: origin.domain()
        Spectrum of Rational Field
        sage: origin.codomain()
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
    """
    pass

