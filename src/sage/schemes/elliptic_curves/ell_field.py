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
from sage.rings.complex_field import is_ComplexField
from sage.rings.real_mpfr import is_RealField
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

    More complicated twists exist in theory for char=2,3 and j=0=1728,
    but are not implemented.
    """

    def quadratic_twist(self, D=None):
        """
        Return the quadratic twist of this curve by ``D``.

        INPUT:

        - ``D`` (default None) the twisting parameter (see below).

        In characteristics other than 2, `D` must be nonzero, and the
        twist is isomorphic to self after adjoining `\sqrt(D)` to the
        base.

        In characteristic 2, `D` is arbitrary, and the twist is
        isomorphic to self after adjoining a root of `x^2+x+D` to the
        base.

        In characteristic 2 when `j=0`, this is not implemented.

        If the base field `F` is finite, `D` need not be specified,
        and the curve returned is the unique curve (up to isomorphism)
        defined over `F` isomorphic to the original curve over the
        quadratic extension of `F` but not over `F` itself.  Over
        infinite fields, an error is raised if `D` is not given.

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
                raise ValueError("twisting parameter D must be specified over infinite fields.")
        else:
            try:
                D=K(D)
            except ValueError:
                raise ValueError("twisting parameter D must be in the base field.")

            if char!=2 and D.is_zero():
                raise ValueError("twisting parameter D must be nonzero when characteristic is not 2")

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
            raise ValueError("Quadratic twist not implemented in char 2 when j=0")

    def two_torsion_rank(self):
        r"""
        Return the dimension of the 2-torsion subgroup of
        `E(K)`.

        This will be 0, 1 or 2.

        EXAMPLES::

            sage: E=EllipticCurve('11a1')
            sage: E.two_torsion_rank()
            0
            sage: K.<alpha>=QQ.extension(E.division_polynomial(2).monic())
            sage: E.base_extend(K).two_torsion_rank()
            1
            sage: E.reduction(53).two_torsion_rank()
            2

        ::

            sage: E = EllipticCurve('14a1')
            sage: E.two_torsion_rank()
            1
            sage: K.<alpha>=QQ.extension(E.division_polynomial(2).monic().factor()[1][0])
            sage: E.base_extend(K).two_torsion_rank()
            2

        ::

            sage: EllipticCurve('15a1').two_torsion_rank()
            2

        """
        f=self.division_polynomial(rings.Integer(2))
        n=len(f.roots())+1
        return rings.Integer(n).ord(rings.Integer(2))


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
            raise ValueError("Quartic twist not defined in chars 2,3")

        if self.j_invariant() !=K(1728):
            raise ValueError("Quartic twist not defined when j!=1728")

        if D.is_zero():
            raise ValueError("quartic twist requires a nonzero argument")

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
            raise ValueError("Sextic twist not defined in chars 2,3")

        if self.j_invariant() !=K(0):
            raise ValueError("Sextic twist not defined when j!=0")

        if D.is_zero():
            raise ValueError("Sextic twist requires a nonzero argument")

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

        If the curves are defined over `\mathbb{Q}`, the output `D` is
        a squarefree integer.

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
            sage: E1.is_quadratic_twist(E1)
            1
            sage: type(E1.is_quadratic_twist(E1)) == type(E1.is_quadratic_twist(E2))   #trac 6574
            True

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
            raise ValueError("arguments are not elliptic curves")
        K = E.base_ring()
        zero = K.zero_element()
        if not K == F.base_ring():
            return zero
        j=E.j_invariant()
        if  j != F.j_invariant():
            return zero

        if E.is_isomorphic(F):
            if K is rings.QQ:
                return rings.ZZ(1)
            return K.one_element()

        char=K.characteristic()

        if char==2:
            raise NotImplementedError("not implemented in characteristic 2")
        elif char==3:
            if j==0:
                raise NotImplementedError("not implemented in characteristic 3 for curves of j-invariant 0")
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
            raise ValueError("arguments are not elliptic curves")
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
            raise NotImplementedError("not implemented in characteristic 2")
        elif char==3:
            raise NotImplementedError("not implemented in characteristic 3")
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
            raise ValueError("arguments are not elliptic curves")
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
            raise NotImplementedError("not implemented in characteristic 2")
        elif char==3:
            raise NotImplementedError("not implemented in characteristic 3")
        else:
            # now char!=2,3:
            D = F.c6()/E.c6()

        if D.is_zero():
            return D

        assert E.sextic_twist(D).is_isomorphic(F)

        return D

    def descend_to(self, K, f=None):
        r"""
        Given a subfield `K` and an elliptic curve self defined over a field `L`,
        this function determines whether there exists an elliptic curve over `K`
        which is isomorphic over `L` to self. If one exists, it finds it.

        INPUT:

        - `K` -- a subfield of the base field of self.
        - `f` -- an embedding of `K` into the base field of self.

        OUTPUT:

        Either an elliptic curve defined over `K` which is isomorphic to self
        or None if no such curve exists.

        .. NOTE::

            This only works over number fields and QQ.

        EXAMPLES::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.descend_to(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Input must be a field.

        ::

            sage: F.<b> = QuadraticField(23)
            sage: G.<a> = F.extension(x^3+5)
            sage: E = EllipticCurve(j=1728*b).change_ring(G)
            sage: E.descend_to(F)
            Elliptic Curve defined by y^2 = x^3 + (8957952*b-206032896)*x + (-247669456896*b+474699792384) over Number Field in b with defining polynomial x^2 - 23

        ::

            sage: L.<a> = NumberField(x^4 - 7)
            sage: K.<b> = NumberField(x^2 - 7)
            sage: E = EllipticCurve([a^6,0])
            sage: E.descend_to(K)
            Elliptic Curve defined by y^2 = x^3 + 1296/49*b*x over Number Field in b with defining polynomial x^2 - 7

        ::

            sage: K.<a> = QuadraticField(17)
            sage: E = EllipticCurve(j = 2*a)
            sage: print E.descend_to(QQ)
            None
        """
        if not K.is_field():
            raise TypeError("Input must be a field.")
        if self.base_field()==K:
            return self
        j = self.j_invariant()
        from sage.rings.all import QQ
        if K == QQ:
            f = QQ.embeddings(self.base_field())[0]
            if j in QQ:
                jbase = QQ(j)
            else:
                return None
        elif f is None:
            embeddings = K.embeddings(self.base_field())
            if len(embeddings) == 0:
                raise TypeError("Input must be a subfield of the base field of the curve.")
            for g in embeddings:
                try:
                    jbase = g.preimage(j)
                    f = g
                    break
                except Exception:
                    pass
            if f is None:
                return None
        else:
            try:
                jbase = f.preimage(j)
            except Exception:
                return None
        E = EllipticCurve(j=jbase)
        E2 = EllipticCurve(self.base_field(), [f(a) for a in E.a_invariants()])
        if jbase==0:
            d = self.is_sextic_twist(E2)
            if d == 1:
                return E
            if d == 0:
                return None
            Etwist = E2.sextic_twist(d)
        elif jbase==1728:
            d = self.is_quartic_twist(E2)
            if d == 1:
                return E
            if d == 0:
                return None
            Etwist = E2.quartic_twist(d)
        else:
            d = self.is_quadratic_twist(E2)
            if d == 1:
                return E
            if d == 0:
                return None
            Etwist = E2.quadratic_twist(d)
        if Etwist.is_isomorphic(self):
            try:
                Eout = EllipticCurve(K, [f.preimage(a) for a in Etwist.a_invariants()])
            except Exception:
                return None
            else:
                return Eout

    def isogeny(self, kernel, codomain=None, degree=None, model=None, check=True):
        r"""
        Returns an elliptic curve isogeny from self.

        The isogeny can be determined in two ways, either by a
        polynomial or a set of torsion points.  The methods used are:

        - Velu's Formulas: Velu's original formulas for computing
          isogenies.  This algorithm is selected by giving as the
          ``kernel`` parameter a point or a list of points which
          generate a finite subgroup.

        - Kohel's Formulas: Kohel's original formulas for computing
          isogenies.  This algorithm is selected by giving as the
          ``kernel`` parameter a polynomial (or a coefficient list
          (little endian)) which will define the kernel of the
          isogeny.

        INPUT:

        - ``E``         - an elliptic curve, the domain of the isogeny to
                          initialize.

        - ``kernel``    - a kernel, either a point in ``E``, a list of points
                          in ``E``, a univariate kernel polynomial or ``None``.
                          If initiating from a domain/codomain, this must be
                          set to None.  Validity of input is *not* fully checked.

        - ``codomain``  - an elliptic curve (default:None).  If ``kernel`` is
                          None, then this must be the codomain of a separable
                          normalized isogeny, furthermore, ``degree`` must be
                          the degree of the isogeny from ``E`` to ``codomain``.
                          If ``kernel`` is not None, then this must be
                          isomorphic to the codomain of the normalized separable
                          isogeny defined by ``kernel``, in this case, the
                          isogeny is post composed with an isomorphism so that
                          this parameter is the codomain.

        - ``degree``    - an integer (default:None). If ``kernel`` is None,
                          then this is the degree of the isogeny from ``E`` to
                          ``codomain``. If ``kernel`` is not None, then this is
                          used to determine whether or not to skip a gcd of the
                          kernel polynomial with the two torsion polynomial of
                          ``E``.

        - ``model``     - a string (default:None).  Only supported variable is
                          "minimal", in which case if``E`` is a curve over the
                          rationals, then the codomain is set to be the unique
                          global minimum model.

        - ``check`` (default: True) does some partial checks that the
                          input is valid (e.g., that the points
                          defined by the kernel polynomial are
                          torsion); however, invalid input can in some
                          cases still pass, since that the points define
                          a group is not checked.

        OUTPUT:

        An isogeny between elliptic curves. This is a morphism of curves.

        EXAMPLES::

            sage: F = GF(2^5, 'alpha'); alpha = F.gen()
            sage: E = EllipticCurve(F, [1,0,1,1,1])
            sage: R.<x> = F[]
            sage: phi = E.isogeny(x+1)
            sage: phi.rational_maps()
            ((x^2 + x + 1)/(x + 1), (x^2*y + x)/(x^2 + 1))

            sage: E = EllipticCurve('11a1')
            sage: P = E.torsion_points()[1]
            sage: E.isogeny(P)
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field

            sage: E = EllipticCurve(GF(19),[1,1])
            sage: P = E(15,3); Q = E(2,12);
            sage: (P.order(), Q.order())
            (7, 3)
            sage: phi = E.isogeny([P,Q]); phi
            Isogeny of degree 21 from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19 to Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19
            sage: phi(E.random_point()) # all points defined over GF(19) are in the kernel
            (0 : 1 : 0)


            # not all polynomials define a finite subgroup trac #6384
            sage: E = EllipticCurve(GF(31),[1,0,0,1,2])
            sage: phi = E.isogeny([14,27,4,1])
            Traceback (most recent call last):
            ...
            ValueError: The polynomial does not define a finite subgroup of the elliptic curve.

        An example in which we construct an invalid morphism, which
        illustrates that the check for correctness of the input is not
        sufficient. (See trac 11578.)::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2-x-1)
            sage: E = EllipticCurve(K, [-13392, -1080432])
            sage: R.<x> = K[]
            sage: phi = E.isogeny( (x-564)*(x - 396/5*a + 348/5) )
            sage: phi.codomain().conductor().norm().factor()
            5^2 * 11^2 * 3271 * 15806939 * 4169267639351
            sage: phi.domain().conductor().norm().factor()
            11^2
        """
        return EllipticCurveIsogeny(self, kernel, codomain, degree, model, check=check)


    def isogeny_codomain(self, kernel, degree=None):
        r"""
        Returns the codomain of the isogeny from self with given
        kernel.

        INPUT:

        - ``kernel`` - Either a list of points in the kernel of the isogeny,
                       or a kernel polynomial (specified as a either a
                       univariate polynomial or a coefficient list.)

        - ``degree`` - an integer, (default:None) optionally specified degree
                       of the kernel.

        OUTPUT:

        An elliptic curve, the codomain of the separable normalized
        isogeny from this kernel

        EXAMPLES::

            sage: E = EllipticCurve('17a1')
            sage: R.<x> = QQ[]
            sage: E2 = E.isogeny_codomain(x - 11/4); E2
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 1461/16*x - 19681/64 over Rational Field

        """
        return isogeny_codomain_from_kernel(self, kernel, degree=None)

    def isogenies_prime_degree(self, l=None, max_l=31):
        """
        Generic code, valid for all fields, for arbitrary prime `l` not equal to the characteristic.

        INPUT:

        - ``l`` -- either None, a prime or a list of primes.
        - ``max_l`` -- a bound on the primes to be tested (ignored unless `l` is None).

        OUTPUT:

        (list) All `l`-isogenies for the given `l` with domain self.

        METHOD:

        Calls the generic function
        ``isogenies_prime_degree()``.  This requires that
        certain operations have been implemented over the base field,
        such as root-finding for univariate polynomials.

        EXAMPLES::

            sage: F = QQbar
            sage: E = EllipticCurve(F, [1,18]); E
            Elliptic Curve defined by y^2 = x^3 + x + 18 over Algebraic Field
            sage: E.isogenies_prime_degree()
            Traceback (most recent call last):
            ...
            NotImplementedError: This code could be implemented for QQbar, but has not been yet.

            sage: F = CC
            sage: E = EllipticCurve(F, [1,18]); E
            Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 18.0000000000000 over Complex Field with 53 bits of precision
            sage: E.isogenies_prime_degree(11)
            Traceback (most recent call last):
            ...
            NotImplementedError: This code could be implemented for general complex fields, but has not been yet.

        Examples over finite fields::

            sage: E = EllipticCurve(GF(next_prime(1000000)), [7,8])
            sage: E.isogenies_prime_degree()
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003, Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003, Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003, Isogeny of degree 17 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 347438*x + 594729 over Finite Field of size 1000003, Isogeny of degree 17 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 674846*x + 7392 over Finite Field of size 1000003, Isogeny of degree 23 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 390065*x + 605596 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(3)
            []
            sage: E.isogenies_prime_degree(5)
            []
            sage: E.isogenies_prime_degree(7)
            []
            sage: E.isogenies_prime_degree(13)
            [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
            Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003]

            sage: E.isogenies_prime_degree([2, 3, 5, 7, 13])
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003, Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003, Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree([2, 4])
            Traceback (most recent call last):
            ...
            ValueError: 4 is not prime.
            sage: E.isogenies_prime_degree(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not prime.
            sage: E.isogenies_prime_degree(11)
            []
            sage: E = EllipticCurve(GF(17),[2,0])
            sage: E.isogenies_prime_degree(3)
            []
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17 to Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 17, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17 to Elliptic Curve defined by y^2 = x^3 + 5*x + 9 over Finite Field of size 17, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17 to Elliptic Curve defined by y^2 = x^3 + 5*x + 8 over Finite Field of size 17]

            sage: E = EllipticCurve(GF(13^4, 'a'),[2,8])
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in a of size 13^4 to Elliptic Curve defined by y^2 = x^3 + 7*x + 4 over Finite Field in a of size 13^4, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in a of size 13^4 to Elliptic Curve defined by y^2 = x^3 + (8*a^3+2*a^2+7*a+5)*x + (12*a^3+3*a^2+4*a+4) over Finite Field in a of size 13^4, Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in a of size 13^4 to Elliptic Curve defined by y^2 = x^3 + (5*a^3+11*a^2+6*a+11)*x + (a^3+10*a^2+9*a) over Finite Field in a of size 13^4]

            sage: E.isogenies_prime_degree(3)
            [Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in a of size 13^4 to Elliptic Curve defined by y^2 = x^3 + 9*x + 11 over Finite Field in a of size 13^4]

        Example to show that separable isogenies of degree equal to the characteristic are now implemented::

            sage: E.isogenies_prime_degree(13)
            [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in a of size 13^4 to Elliptic Curve defined by y^2 = x^3 + 6*x + 5 over Finite Field in a of size 13^4]

        Examples over number fields (other than QQ)::

            sage: QQroot2.<e> = NumberField(x^2-2)
            sage: E = EllipticCurve(QQroot2, j=8000)
            sage: E.isogenies_prime_degree()
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 = x^3 + (-602112000)*x + 5035261952000 over Number Field in e with defining polynomial x^2 - 2,
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 = x^3 + (903168000*e-1053696000)*x + (14161674240000*e-23288086528000) over Number Field in e with defining polynomial x^2 - 2,
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 = x^3 + (-903168000*e-1053696000)*x + (-14161674240000*e-23288086528000) over Number Field in e with defining polynomial x^2 - 2]

            sage: E = EllipticCurve(QQroot2, [1,0,1,4, -6]); E
            Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-36)*x + (-70) over Number Field in e with defining polynomial x^2 - 2]
            sage: E.isogenies_prime_degree(3)
            [Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-128/3)*x + 5662/27 over Number Field in e with defining polynomial x^2 - 2, Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-171)*x + (-874) over Number Field in e with defining polynomial x^2 - 2]
        """
        F = self.base_ring()
        if is_RealField(F):
            raise NotImplementedError("This code could be implemented for general real fields, but has not been yet.")
        if is_ComplexField(F):
            raise NotImplementedError("This code could be implemented for general complex fields, but has not been yet.")
        if F == rings.QQbar:
            raise NotImplementedError("This code could be implemented for QQbar, but has not been yet.")

        from isogeny_small_degree import isogenies_prime_degree
        if l is None:
            from sage.rings.all import prime_range
            l = prime_range(max_l+1)

        if not isinstance(l, list):
            try:
                l = rings.ZZ(l)
            except TypeError:
                raise ValueError("%s is not prime."%l)
            if l.is_prime():
                return isogenies_prime_degree(self, l)
            else:
                raise ValueError("%s is not prime."%l)

        L = list(set(l))
        try:
            L = [rings.ZZ(l) for l in L]
        except TypeError:
            raise ValueError("%s is not a list of primes."%l)

        L.sort()
        return sum([isogenies_prime_degree(self,l) for l in L],[])

    def is_isogenous(self, other, field=None):
        """
        Returns whether or not self is isogenous to other.

        INPUT:

        - ``other`` -- another elliptic curve.

        - ``field`` (default None) -- Currently not implemented. A
          field containing the base fields of the two elliptic curves
          onto which the two curves may be extended to test if they
          are isogenous over this field. By default is_isogenous will
          not try to find this field unless one of the curves can be
          be extended into the base field of the other, in which case
          it will test over the larger base field.

        OUTPUT:

        (bool) True if there is an isogeny from curve ``self`` to
        curve ``other`` defined over ``field``.

        METHOD:

        Over general fields this is only implemented in trivial cases.

        EXAMPLES::

            sage: E1 = EllipticCurve(CC, [1,18]); E1
            Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 18.0000000000000 over Complex Field with 53 bits of precision
            sage: E2 = EllipticCurve(CC, [2,7]); E2
            Elliptic Curve defined by y^2 = x^3 + 2.00000000000000*x + 7.00000000000000 over Complex Field with 53 bits of precision
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented for isomorphic curves over general fields.

            sage: E1 = EllipticCurve(Frac(PolynomialRing(ZZ,'t')), [2,19]); E1
            Elliptic Curve defined by y^2 = x^3 + 2*x + 19 over Fraction Field of Univariate Polynomial Ring in t over Integer Ring
            sage: E2 = EllipticCurve(CC, [23,4]); E2
            Elliptic Curve defined by y^2 = x^3 + 23.0000000000000*x + 4.00000000000000 over Complex Field with 53 bits of precision
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented for isomorphic curves over general fields.
        """
        from ell_generic import is_EllipticCurve
        if not is_EllipticCurve(other):
            raise ValueError("Second argument is not an Elliptic Curve.")
        if self.is_isomorphic(other):
            return True
        else:
            raise NotImplementedError("Only implemented for isomorphic curves over general fields.")

    def weierstrass_p(self, prec=20, algorithm=None):
        r"""
        Computes the Weierstrass `\wp`-function of the elliptic curve.

        INPUT:

        - ``mprec`` - precision

        - ``algorithm`` - string (default:``None``) an algorithm identifier
                      indicating using the ``pari``, ``fast`` or ``quadratic``
                      algorithm. If the algorithm is ``None``, then this
                      function determines the best algorithm to use.

        OUTPUT:

        a Laurent series in one variable `z` with coefficients in the
        base field `k` of `E`.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: E.weierstrass_p(prec=10)
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + O(z^10)
            sage: E.weierstrass_p(prec=8)
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + O(z^8)
            sage: Esh = E.short_weierstrass_model()
            sage: Esh.weierstrass_p(prec=8)
            z^-2 + 13392/5*z^2 + 1080432/7*z^4 + 59781888/25*z^6 + O(z^8)
            sage: E.weierstrass_p(prec=20, algorithm='fast')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
            sage: E.weierstrass_p(prec=20, algorithm='pari')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
            sage: E.weierstrass_p(prec=20, algorithm='quadratic')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
        """
        from ell_wp import weierstrass_p
        return weierstrass_p(self, prec=prec, algorithm=algorithm)

    def hasse_invariant(self):
        r"""
        Returns the Hasse invariant of this elliptic curve.

        OUTPUT:

        The Hasse invariant of this elliptic curve, as an element of
        the base field.  This is only defined over fields of positive
        characteristic, and is an element of the field which is zero
        if and only if the curve is supersingular.  Over a field of
        characteristic zero, where the Hasse invariant is undefined,
        a ``ValueError`` is returned.

        EXAMPLES::

            sage: E = EllipticCurve([Mod(1,2),Mod(1,2),0,0,Mod(1,2)])
            sage: E.hasse_invariant()
            1
            sage: E = EllipticCurve([0,0,Mod(1,3),Mod(1,3),Mod(1,3)])
            sage: E.hasse_invariant()
            0
            sage: E = EllipticCurve([0,0,Mod(1,5),0,Mod(2,5)])
            sage: E.hasse_invariant()
            0
            sage: E = EllipticCurve([0,0,Mod(1,5),Mod(1,5),Mod(2,5)])
            sage: E.hasse_invariant()
            2

        Some examples over larger fields::

            sage: EllipticCurve(GF(101),[0,0,0,0,1]).hasse_invariant()
            0
            sage: EllipticCurve(GF(101),[0,0,0,1,1]).hasse_invariant()
            98
            sage: EllipticCurve(GF(103),[0,0,0,0,1]).hasse_invariant()
            20
            sage: EllipticCurve(GF(103),[0,0,0,1,1]).hasse_invariant()
            17
            sage: F.<a> = GF(107^2)
            sage: EllipticCurve(F,[0,0,0,a,1]).hasse_invariant()
            62*a + 75
            sage: EllipticCurve(F,[0,0,0,0,a]).hasse_invariant()
            0

        Over fields of characteristic zero, the Hasse invariant is
        undefined::

            sage: E = EllipticCurve([0,0,0,0,1])
            sage: E.hasse_invariant()
            Traceback (most recent call last):
            ...
            ValueError: Hasse invariant only defined in positive characteristic
        """
        k = self.base_field()
        p = k.characteristic()
        if p == 0:
            raise ValueError('Hasse invariant only defined in positive characteristic')
        elif p == 2:
            return self.a1()
        elif p == 3:
            return self.b2()
        elif p == 5:
            return self.c4()
        elif p == 7:
            return -self.c6()
        else:
            R = k['x']
            x = R.gen()
            E = self.short_weierstrass_model()
            f=(x**3+E.a4()*x+E.a6())**((p-1)//2)
            return f.coeffs()[p-1]
