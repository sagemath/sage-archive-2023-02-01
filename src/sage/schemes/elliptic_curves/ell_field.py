r"""
Elliptic curves over a general field

This module defines the class ``EllipticCurve_field``, based on
``EllipticCurve_generic``, for elliptic curves over general fields.

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import ell_generic
import sage.rings.all as rings
from constructor import EllipticCurve

from ell_curve_isogeny import EllipticCurveIsogeny, isogeny_codomain_from_kernel

class EllipticCurve_field(ell_generic.EllipticCurve_generic):

    base_field = ell_generic.EllipticCurve_generic.base_ring

    # Twists: rewritten by John Cremona as follows:
    #
    # Quadratic twist allowed except when char=2, j=0
    # Quartic twist allowed only if j=1728!=0 (so char!=2,3)
    # Sextic  twist allowed only if j=0!=1728 (so char!=2,3)
    #
    # More complicated twists exist in theory for char=2,3 and
    # j=0=1728, but I have never worked them out or seen them used!
    #

    r"""
    Twists: rewritten by John Cremona as follows:

    The following twists are implemented:

    - Quadratic twist:  except when char=2 and `j=0`.
    - Quartic twist: only if `j=1728\not=0` (so not if char=2,3).
    - Sextic  twist: only if `j=0\not=1728` (so not if char=2,3).

    More complicated twists exist in theory for char=2,3 and
    j=0=1728, but are not implemented.
    """

    def quadratic_twist(self, D=None):
        """
        Return the quadratic twist of this curve by ``D``.

        INPUT:

        - ``D`` (default None) the twisting parameter (see below).

        In characteristics other than 2, `D` must be nonzero, and the twist is
        isomorphic to self after adjoining `\sqrt(D)` to the base.

        In characteristic 2, `D` is arbitrary, and the twist is isomorphic
        to self after adjoining a root of `x^2+x+D` to the base.

        In characteristic 2 when `j=0`, this is not implemented.

        If the base field `F` is finite, `D` need not be specified, and
        the curve returned is the unique curve (up to isomorphism)
        defined over `F` isomorphic to the original curve over the
        quadratic extension of `F` but not over `F` itself.  Over infinte
        fields, an error is raised if `D` is not given.

        EXAMPLES::

            sage: E = EllipticCurve([GF(1103)(1), 0, 0, 107, 340]); E
            Elliptic Curve defined by y^2 + x*y  = x^3 + 107*x + 340 over Finite Field of size 1103
            sage: F=E.quadratic_twist(-1); F
            Elliptic Curve defined by y^2  = x^3 + 1102*x^2 + 609*x + 300 over Finite Field of size 1103
            sage: E.is_isomorphic(F)
            False
            sage: E.is_isomorphic(F,GF(1103^2,'a'))
            True

        A characteristic 2 example::

            sage: E=EllipticCurve(GF(2),[1,0,1,1,1])
            sage: E1=E.quadratic_twist(1)
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1,GF(4,'a'))
            True

        Over finite fields, the twisting parameter may be omitted::

            sage: k.<a> = GF(2^10)
            sage: E = EllipticCurve(k,[a^2,a,1,a+1,1])
            sage: Et = E.quadratic_twist()
            sage: Et # random (only determined up to isomorphism)
            Elliptic Curve defined by y^2 + x*y  = x^3 + (a^7+a^4+a^3+a^2+a+1)*x^2 + (a^8+a^6+a^4+1) over Finite Field in a of size 2^10
            sage: E.is_isomorphic(Et)
            False
            sage: E.j_invariant()==Et.j_invariant()
            True

            sage: p=next_prime(10^10)
            sage: k = GF(p)
            sage: E = EllipticCurve(k,[1,2,3,4,5])
            sage: Et = E.quadratic_twist()
            sage: Et # random (only determined up to isomorphism)
            Elliptic Curve defined by y^2  = x^3 + 7860088097*x^2 + 9495240877*x + 3048660957 over Finite Field of size 10000000019
            sage: E.is_isomorphic(Et)
            False
            sage: k2 = GF(p^2,'a')
            sage: E.change_ring(k2).is_isomorphic(Et.change_ring(k2))
            True
        """
        K=self.base_ring()
        char=K.characteristic()

        if D is None:
            if K.is_finite():
                x = rings.polygen(K)
                if char==2:
                    # We find D such that x^2+x+D is irreducible. If the
                    # degree is odd we can take D=1; otherwise it suffices to
                    # consider odd powers of a generator.
                    D = K(1)
                    if K.degree()%2==0:
                        D = K.gen()
                        a = D**2
                        while len((x**2+x+D).roots())>0:
                            D *= a
                else:
                    # We could take a multiplicative generator but
                    # that might be expensive to compute; otherwise
                    # half the elements will do
                    D = K.random_element()
                    while len((x**2-D).roots())>0:
                        D = K.random_element()
            else:
                raise ValueError, "twisting parameter D must be specified over infinite fields."
        else:
            try:
                D=K(D)
            except ValueError:
                raise ValueError, "twisting parameter D must be in the base field."

            if char!=2 and D.is_zero():
                raise ValueError, "twisting parameter D must be nonzero when characteristic is not 2"

        if char!=2:
            b2,b4,b6,b8=self.b_invariants()
            # E is isomorphic to  [0,b2,0,8*b4,16*b6]
            return EllipticCurve(K,[0,b2*D,0,8*b4*D**2,16*b6*D**3])

        # now char==2
        if self.j_invariant() !=0: # iff a1!=0
            a1,a2,a3,a4,a6=self.ainvs()
            E0=self.change_weierstrass_model(a1,a3/a1,0,(a1**2*a4+a3**2)/a1**3)
            # which has the form = [1,A2,0,0,A6]
            assert E0.a1()==K(1)
            assert E0.a3()==K(0)
            assert E0.a4()==K(0)
            return EllipticCurve(K,[1,E0.a2()+D,0,0,E0.a6()])
        else:
            raise ValueError, "Quadratic twist not implemented in char 2 when j=0"

    def quartic_twist(self, D):
        r"""
        Return the quartic twist of this curve by `D`.

        INPUT:

        - ``D`` (must be nonzero) -- the twisting parameter..

        .. note::

           The characteristic must not be 2 or 3, and the `j`-invariant must be 1728.

        EXAMPLES::

            sage: E=EllipticCurve_from_j(GF(13)(1728)); E
            Elliptic Curve defined by y^2  = x^3 + x over Finite Field of size 13
            sage: E1=E.quartic_twist(2); E1
            Elliptic Curve defined by y^2  = x^3 + 5*x over Finite Field of size 13
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1,GF(13^2,'a'))
            False
            sage: E.is_isomorphic(E1,GF(13^4,'a'))
            True
        """
        K=self.base_ring()
        char=K.characteristic()
        D=K(D)

        if char==2 or char==3:
            raise ValueError, "Quartic twist not defined in chars 2,3"

        if self.j_invariant() !=K(1728):
            raise ValueError, "Quartic twist not defined when j!=1728"

        if D.is_zero():
            raise ValueError, "quartic twist requires a nonzero argument"

        c4,c6=self.c_invariants()
        # E is isomorphic to  [0,0,0,-27*c4,0]
        assert c6==0
        return EllipticCurve(K,[0,0,0,-27*c4*D,0])

    def sextic_twist(self, D):
        r"""
        Return the quartic twist of this curve by `D`.

        INPUT:

        - ``D`` (must be nonzero) -- the twisting parameter..

        .. note::

           The characteristic must not be 2 or 3, and the `j`-invariant must be 0.

        EXAMPLES::

            sage: E=EllipticCurve_from_j(GF(13)(0)); E
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 13
            sage: E1=E.sextic_twist(2); E1
            Elliptic Curve defined by y^2 = x^3 + 11 over Finite Field of size 13
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1,GF(13^2,'a'))
            False
            sage: E.is_isomorphic(E1,GF(13^4,'a'))
            False
            sage: E.is_isomorphic(E1,GF(13^6,'a'))
            True
        """
        K=self.base_ring()
        char=K.characteristic()
        D=K(D)

        if char==2 or char==3:
            raise ValueError, "Sextic twist not defined in chars 2,3"

        if self.j_invariant() !=K(0):
            raise ValueError, "Sextic twist not defined when j!=0"

        if D.is_zero():
            raise ValueError, "Sextic twist requires a nonzero argument"

        c4,c6=self.c_invariants()
        # E is isomorphic to  [0,0,0,0,-54*c6]
        assert c4==0
        return EllipticCurve(K,[0,0,0,0,-54*c6*D])

    def is_quadratic_twist(self, other):
        r"""
        Determine whether this curve is a quadratic twist of another.

        INPUT:

        - ``other`` -- an elliptic curves with the same base field as self.

        OUTPUT:

        Either 0, if the curves are not quadratic twists, or `D` if
        ``other`` is ``self.quadratic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        .. note::

           Not fully implemented in characteristic 2, or in
           characteristic 3 when both `j`-invariants are 0.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Et = E.quadratic_twist(-24)
            sage: E.is_quadratic_twist(Et)
            -6

            sage: E1=EllipticCurve([0,0,1,0,0])
            sage: E1.j_invariant()
            0
            sage: E2=EllipticCurve([0,0,0,0,2])
            sage: E1.is_quadratic_twist(E2)
            2

        ::

            sage: E1=EllipticCurve([0,0,0,1,0])
            sage: E1.j_invariant()
            1728
            sage: E2=EllipticCurve([0,0,0,2,0])
            sage: E1.is_quadratic_twist(E2)
            0
            sage: E2=EllipticCurve([0,0,0,25,0])
            sage: E1.is_quadratic_twist(E2)
            5

        ::

            sage: F = GF(101)
            sage: E1 = EllipticCurve(F,[4,7])
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2); D!=0
            True
            sage: F = GF(101)
            sage: E1 = EllipticCurve(F,[4,7])
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2)
            sage: E1.quadratic_twist(D).is_isomorphic(E2)
            True
            sage: E1.is_isomorphic(E2)
            False
            sage: F2 = GF(101^2,'a')
            sage: E1.change_ring(F2).is_isomorphic(E2.change_ring(F2))
            True

        A characteristic 3 example::

            sage: F = GF(3^5,'a')
            sage: E1 = EllipticCurve_from_j(F(1))
            sage: E2 = E1.quadratic_twist(-1)
            sage: D = E1.is_quadratic_twist(E2); D!=0
            True
            sage: E1.quadratic_twist(D).is_isomorphic(E2)
            True

        ::

            sage: E1 = EllipticCurve_from_j(F(0))
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2); D
            1
            sage: E1.is_isomorphic(E2)
            True

        """
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
            raise ValueError, "arguments are not elliptic curves"
        K = E.base_ring()
        zero = K.zero_element()
        if not K == F.base_ring():
            return zero
        j=E.j_invariant()
        if  j != F.j_invariant():
            return zero

        if E.is_isomorphic(F):
            return K.one_element()

        char=K.characteristic()

        if char==2:
            raise NotImplementedError, "not implemented in characteristic 2"
        elif char==3:
            if j==0:
                raise NotImplementedError, "not implemented in characteristic 3 for curves of j-invariant 0"
            D = E.b2()/F.b2()

        else:
            # now char!=2,3:
            c4E,c6E = E.c_invariants()
            c4F,c6F = F.c_invariants()

            if j==0:
                um = c6E/c6F
                x=rings.polygen(K)
                ulist=(x**3-um).roots(multiplicities=False)
                if len(ulist)==0:
                    D = zero
                else:
                    D = ulist[0]
            elif j==1728:
                um=c4E/c4F
                x=rings.polygen(K)
                ulist=(x**2-um).roots(multiplicities=False)
                if len(ulist)==0:
                    D = zero
                else:
                    D = ulist[0]
            else:
                D = (c6E*c4F)/(c6F*c4E)

        # Normalization of output:

        if D.is_zero():
            return D

        if K is rings.QQ:
            D = D.squarefree_part()

        assert E.quadratic_twist(D).is_isomorphic(F)

        return D

    def is_quartic_twist(self, other):
        r"""
        Determine whether this curve is a quartic twist of another.

        INPUT:

        - ``other`` -- an elliptic curves with the same base field as self.

        OUTPUT:

        Either 0, if the curves are not quartic twists, or `D` if
        ``other`` is ``self.quartic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        .. note::

           Not fully implemented in characteristics 2 or 3.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(1728))
            sage: E1 = E.quartic_twist(2)
            sage: D = E.is_quartic_twist(E1); D!=0
            True
            sage: E.quartic_twist(D).is_isomorphic(E1)
            True

        ::

            sage: E = EllipticCurve_from_j(1728)
            sage: E1 = E.quartic_twist(12345)
            sage: D = E.is_quartic_twist(E1); D
            15999120
            sage: (D/12345).is_perfect_power(4)
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
            raise ValueError, "arguments are not elliptic curves"
        K = E.base_ring()
        zero = K.zero_element()
        if not K == F.base_ring():
            return zero
        j=E.j_invariant()
        if  j != F.j_invariant() or j!=K(1728):
            return zero

        if E.is_isomorphic(F):
            return K.one_element()

        char=K.characteristic()

        if char==2:
            raise NotImplementedError, "not implemented in characteristic 2"
        elif char==3:
            raise NotImplementedError, "not implemented in characteristic 3"
        else:
            # now char!=2,3:
            D = F.c4()/E.c4()

        if D.is_zero():
            return D

        assert E.quartic_twist(D).is_isomorphic(F)

        return D

    def is_sextic_twist(self, other):
        r"""
        Determine whether this curve is a sextic twist of another.

        INPUT:

        - ``other`` -- an elliptic curves with the same base field as self.

        OUTPUT:

        Either 0, if the curves are not sextic twists, or `D` if
        ``other`` is ``self.sextic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        .. note::

           Not fully implemented in characteristics 2 or 3.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(0))
            sage: E1 = E.sextic_twist(2)
            sage: D = E.is_sextic_twist(E1); D!=0
            True
            sage: E.sextic_twist(D).is_isomorphic(E1)
            True

        ::

            sage: E = EllipticCurve_from_j(0)
            sage: E1 = E.sextic_twist(12345)
            sage: D = E.is_sextic_twist(E1); D
            575968320
            sage: (D/12345).is_perfect_power(6)
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
            raise ValueError, "arguments are not elliptic curves"
        K = E.base_ring()
        zero = K.zero_element()
        if not K == F.base_ring():
            return zero
        j=E.j_invariant()
        if  j != F.j_invariant() or not j.is_zero():
            return zero

        if E.is_isomorphic(F):
            return K.one_element()

        char=K.characteristic()

        if char==2:
            raise NotImplementedError, "not implemented in characteristic 2"
        elif char==3:
            raise NotImplementedError, "not implemented in characteristic 3"
        else:
            # now char!=2,3:
            D = F.c6()/E.c6()

        if D.is_zero():
            return D

        assert E.sextic_twist(D).is_isomorphic(F)

        return D


    def isogeny(self, kernel, codomain=None, degree=None, model=None, algorithm=None):
        r"""
        Returns an elliptic curve isogeny from self.
        This redirects to the EllipticCurveIsogeny constructor.
        See the EllipticCurveIsogeny class for more information.

        EXAMPLES:

            sage: F = GF(2^5, 'alpha'); alpha = F.gen()
            sage: E = EllipticCurve(F, [1,0,1,1,1])
            sage: R.<x> = F[]
            sage: phi = E.isogeny(x+1)
            sage: phi.rational_maps()
            ((x^2 + x + 1)/(x + 1), (x^2*y + x)/(x^2 + 1))

        """
        return EllipticCurveIsogeny(self, kernel, codomain, degree, model, algorithm)


    def isogeny_codomain(self, kernel, degree=None, algorithm=None):
        r"""
        Returns the codomain of the isogeny from self with givien kernel.
        This redirects to the isogeny_codomain_from_kernel function,
        see that function for more information.

        EXAMPLES:

            sage: E = EllipticCurve('17a1')
            sage: R.<x> = QQ[]
            sage: E2 = E.isogeny_codomain(x - 11/4); E2
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 1461/16*x - 19681/64 over Rational Field

        """
        return isogeny_codomain_from_kernel(self, kernel, degree=None, algorithm=None)
