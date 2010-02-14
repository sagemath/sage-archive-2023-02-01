# -*- coding: utf-8 -*-
r"""
Galois representations attached to elliptic curves

If `E` is an elliptic curve over a global field `K`, currently only implemented over `\QQ`.

Given an elliptic curve `E` over a number field `K`
and a rational prime number `p`, the `p^n`-torsion
`E[p^n]` points of `E` is a representation of the
absolute Galois group `G_K` of `K`. As `n` varies
we obtain the Tate module `T_p E` which is a
a representation of `G_K` on a free `\ZZ_p`-module
of rank `2`. As `p` varies the representations
are compatible.

Currently sage can decide whether the Galois module
`E[p]` is reducible, i.e. if `E` admits an isogeny
of degree `p`, and whether the image of
the representation on `E[p]` is surjective onto
`Aut(E[p]) = GL_2(\mathbb{F}_p)`.

EXAMPLES::

    sage: E = EllipticCurve('196a1')
    sage: rho = E.galois_representation()
    sage: rho.is_irreducible(7)
    True
    sage: rho.is_reducible(3)
    True
    sage: rho.is_irreducible(2)
    True
    sage: rho.is_surjective(2)
    False
    sage: rho.is_surjective(3)
    False
    sage: rho.is_surjective(5)
    True
    sage: rho.reducible_primes()
    [3]
    sage: rho.non_surjective()
    [2, 3]

For semi-stable curve it is known that the representation is
surjective if and only if it is irreducible::

    sage: E = EllipticCurve('11a1')
    sage: rho = E.galois_representation()
    sage: rho.non_surjective()
    [5]
    sage: rho.reducible_primes()
    [5]


For cm curves it is not true that there are only finitely many primes for which the
Galois representation mod p is surjective onto `GL_2(\mathbb{F}_p)`::

    sage: E = EllipticCurve('27a1')
    sage: rho = E.galois_representation()
    sage: rho.non_surjective()
    [0]
    sage: rho.reducible_primes()
    [3]
    sage: E.has_cm()
    True

REFERENCES:

.. [Se1] Jean-Pierre Serre,
    Propriétés galoisiennes des points d'ordre fini
    des courbes elliptiques.
    Invent. Math.  15  (1972), no. 4, 259--331.
.. [Se2] Jean-Pierre Serre,
    Sur les représentations modulaires de degré
    2 de Gal`(\overline\QQ/\QQ)`.
    Duke Math. J. 54 (1987), no. 1, 179--230.
.. [Co] Alina Carmen Cojocaru,
    On the surjectivity of the Galois representations
    associated to non-CM elliptic curves.
    With an appendix by Ernst Kani.
    Canad. Math. Bull. 48 (2005), no. 1, 16--31.


AUTHORS:

- chris wuthrich (02/10) - moved from ell_rational_field.py.

"""

######################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
######################################################################

from sage.structure.sage_object import SageObject
import sage.rings.arith as arith
import sage.misc.misc as misc
import sage.rings.all as rings
from math import sqrt
from sage.libs.pari.all import pari, PariError

class GaloisRepresentation(SageObject):
    r"""
    The compatible family of Galois representation
    attached to an elliptic curve over a number field.

    Given an elliptic curve `E` over a number field `K`
    and a rational prime number `p`, the `p^n`-torsion
    `E[p^n]` points of `E` is a representation of the
    absolute Galois group `G_K` of `K`. As `n` varies
    we obtain the Tate module `T_p E` which is a
    a representation of `G_K` on a free `\ZZ_p`-module
    of rank `2`. As `p` varies the representations
    are compatible.

    EXAMPLES::

        sage: rho = EllipticCurve('11a1').galois_representation()
        sage: rho
        Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

    """

    def __init__(self, E):
        r"""

        see ``GaloisRepresentation`` for documentation

        EXAMPLES::

            sage: rho = EllipticCurve('11a1').galois_representation()
            sage: loads(rho.dumps()) == rho
            True

        """
        self.E = E

    def __repr__(self):
        r"""
        string representation of the class

        EXAMPLES::

            sage: rho = EllipticCurve([0,1]).galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field

        """
        return "Compatible family of Galois representations associated to the " + repr(self.E)

    def __eq__(self,other):
        r"""
        Compares two Galois representations.
        We define tho compatible families of representations
        attached to elliptic curves to be isomorphic if the curves are equal

        EXAMPLES::

            sage: rho = EllipticCurve('11a1').galois_representation()
            sage: rho2 = EllipticCurve('11a2').galois_representation()
            sage: rho == rho
            True
            sage: rho == rho2
            False
            sage: rho == 34
            False

        """
        # if rho_E = rho_E' then the L-functions agree,
        # so E and E' are isogenous
        # except for p=2
        # the mod p representations will be different
        # for p dividing the degree of the isogeny
        # anyway, there should not be a _compatible_
        # isomorphism between rho and rho' unless E
        # is isomorphic to E'
        # Note that rho can not depend on the Weierstrass model
        if not type(self) == type(other):
            return False
        return self.E.is_isomorphic(other.E)


    def elliptic_curve(self):
        r"""
        The elliptic curve associated to this representation.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: rho = E.galois_representation()
            sage: rho.elliptic_curve() == E
            True

        """
        return self.E

    ##########################################################
    # Galois Representations
    ##########################################################

    def is_reducible(self, p):
        r"""
        Return True if the mod-p representation is
        reducible. This is equivalent to the existence of an
        isogeny of degree `p` from the elliptic curve.

        INPUT:

        -  ``p`` - a prime number

        .. note::

           The answer is cached.

        EXAMPLES::

            sage: rho = EllipticCurve('121a').galois_representation()
            sage: rho.is_reducible(7)
            False
            sage: rho.is_reducible(11)
            True
            sage: EllipticCurve('11a').galois_representation().is_reducible(5)
            True
            sage: rho = EllipticCurve('11a2').galois_representation()
            sage: rho.is_reducible(5)
            True
            sage: EllipticCurve('11a2').torsion_order()
            1

        """
        try:
            return self.__is_reducible[p]
        except AttributeError:
            self.__is_reducible = {}
        except KeyError:
            pass

        if not arith.is_prime(p):
            raise ValueError, 'p (=%s) must be prime'%p
        # we do is_surjective first, since this is
        # much easier than computing isogeny_class
        t = self.is_surjective(p)
        if t == True:
            self.__is_reducible[p] = False
            return False  # definitely not reducible
        isogeny_matrix = self.E.isogeny_class()[ 1 ]
        v = isogeny_matrix.row(0) # first row
        for a in v:
            if a != 0 and a % p == 0:
                self.__is_reducible[p] = True
                return True
        self.__is_reducible[p] = False
        return False

    def is_irreducible(self, p):
        r"""
        Return True if the mod p representation is irreducible.

        EXAMPLES::

            sage: rho = EllipticCurve('37b').galois_representation()
            sage: rho.is_irreducible(2)
            True
            sage: rho.is_irreducible(3)
            False
            sage: rho.is_reducible(2)
            False
            sage: rho.is_reducible(3)
            True
        """
        return not self.is_reducible(p)

    def is_surjective(self, p, A=1000):
        r"""
        Return True if the mod-p representation is
        surjective onto `Aut(E[p]) = GL_2(\mathbb{F}_p)`.

        False if it is not, or None if we were unable to
        determine whether it is or not.

        .. note::

           The answer is cached.

        INPUT:


        -  ``p`` - int (a prime number)

        -  ``A`` - int (a bound on the number of a_p to use)


        OUTPUT:

        boolean. True if the mod-p representation is surjective
          and False if (probably) not

        EXAMPLES::

            sage: rho = EllipticCurve('37b').galois_representation()
            sage: rho.is_surjective(2)
            True
            sage: rho.is_surjective(3)
            False

        REMARKS:

        1. If `p = 5` then the mod-p representation is
           surjective if and only if the p-adic representation is
           surjective. When `p = 2, 3` there are
           counterexamples. See a very recent paper of Elkies for more
           details when `p=3`.

        2. When p = 3 this function always gives the correct result
           irregardless of A, since it explicitly determines the
           p-division polynomial.
        """
        if not arith.is_prime(p):
            raise TypeError, "p (=%s) must be prime."%p
        A = int(A)
        key = (p, A)
        try:
            return self.__is_surjective[key]
        except KeyError:
            pass
        except AttributeError:
            self.__is_surjective = {}

        ans = self._is_surjective(p, A)
        self.__is_surjective[key] = ans
        return ans

    def _is_surjective(self, p, A):
        r"""
        helper function for ``is_surjective``

        EXAMPLES::

            sage: rho = EllipticCurve('37b').galois_representation()
            sage: rho._is_surjective(7,100)
            True

        """
        T = self.E.torsion_subgroup().order()
        if T % p == 0:
            return False   #, "%s-torsion"%p

        if p == 2:
            # E is isomorphic to  [0,b2,0,8*b4,16*b6]
            b2,b4,b6,b8=self.E.b_invariants()
            R = rings.PolynomialRing(self.E.base_ring(), 'x')
            x = R.gen()
            f = x**3 + b2*x**2 + 8*b4*x + 16*b6
            if not f.is_irreducible():
                return False    #, '2-torsion'
            if arith.is_square(f.discriminant()):
                return False    #, "A3"
            return True  #, None

        if p == 3:
            # Algorithm: Let f be the 3-division polynomial, which is
            # a polynomial of degree 4.  Then I claim that this
            # polynomial has Galois group S_4 if and only if the
            # representation rhobar_{E,3} is surjective.  If the group
            # is S_4, then S_4 is a quotient of the image of
            # rhobar_{E,3}.  Since S_4 has order 24 and GL_2(F_3)
            # has order 48, the only possibility we have to consider
            # is that the image of rhobar is isomorphic to S_4.
            # But this is not the case because S_4 is not a subgroup
            # of GL_2(F_3).    If it were, it would be normal, since
            # it would have index 2.  But there is a *unique* normal
            # subgroup of GL_2(F_3) of index 2, namely SL_2(F_3),
            # and SL_2(F_3) is not isomorphic to S_4 (S_4 has a normal
            # subgroup of index 2 and SL_2(F_3) does not.)
            # (What's a simple way to see that SL_2(F_3) is the
            # unique index-2 normal subgroup?  I didn't see an obvious
            # reason, so just used the NormalSubgroups command in MAGMA
            # and it output exactly one of index 2.)

            # Here's Noam Elkies proof for the other direction:

            #> Let E be an elliptic curve over Q.  Is the mod-3
            #> representation E[3]  surjective if and only if the
            #> (degree 4) division polynomial has Galois group S_4?  I
            #> can see why the group being S_4 implies the
            #> representation is surjective, but the converse is not
            #> clear to me.
            # I would have thought that this is the easier part: to
            # say that E[3] is surjective is to say the 3-torsion
            # field Q(E[3]) has Galois group GL_2(Z/3) over Q.  Let
            # E[3]+ be the subfield fixed by the element -1 of
            # GL_2(Z/3).  Then E[3] has Galois group PGL_2(Z/3), which
            # is identified with S_4 by its action on the four
            # 3-element subgroups of E[3].  Each such subgroup is in
            # turn determined by the x-coordinate shared by its two
            # nonzero points.  So, if E[3] is surjective then any
            # permutation of those x-coordinates is realized by some
            # element of Gal(E[3]+/Q).  Thus the Galois group of the
            # division polynomial (whose roots are those
            # x-coordinates) maps surjectively to S_4, which means it
            # equals S_4.


            f = self.E.division_polynomial(3)
            if not f.is_irreducible():
                return False   #, "reducible_3-divpoly"
            n = pari(f).polgalois()[0]
            if n == 24:
                return True   #, None
            else:
                return False   #, "3-divpoly_galgroup_order_%s"%n

        if self.E.has_cm():
            return False   #, "CM"
        an = self.E.anlist(A)
        ell = 0
        Np = self.E.conductor() * p
        signs = []
        while True:
            ell = arith.next_prime(ell)
            if ell >= A: break
            if Np % ell != 0:
                a_ell = an[int(ell)]
                if a_ell % p != 0:
                    s = arith.kronecker(a_ell**2 - 4*ell, p)
                    #print ell, s
                    if s == 0: continue
                    if not (s in signs):
                        signs.append(s)
                        if len(signs) == 2:
                            return True  #, None

        # could do something further here...
        return False    #, signs

    #def is_semistable(self):

    def reducible_primes(self):
        r"""
        Returns a list of the primes `p` such that the mod
        `p` representation is reducible. For
        all other primes the representation is irreducible.

        EXAMPLES::

            sage: rho = EllipticCurve('225a').galois_representation()
            sage: rho.reducible_primes()
            [3]
        """
        try:
            return self.__reducible_primes
        except AttributeError:
            pass
        C, I = self.E.isogeny_class(algorithm='sage')
        X = set(I.list())
        R = [p for p in X if arith.is_prime(p)]
        self.__reducible_primes = R
        return R

    def non_surjective(self, A=1000):
        r"""
        Returns a list of primes p such that the mod-p representation
        *might* not be surjective (this list
        usually contains 2, because of shortcomings of the algorithm). If `p`
        is not in the returned list, then the mod-p representation
        is provably surjective.

        By a theorem of Serre, there are only finitely
        many primes in this list, except when the curve has
        complex multiplication.

        If the curve has CM, we simply return the
        sequence [0] and do no further computation.

        INPUT:

        -  ``A`` - an integer

        OUTPUT:

        -  ``list`` - if curve has CM, returns [0].
           Otherwise, returns a list of primes where mod-p representation is
           very likely not surjective. At any prime not in this list, the
           representation is definitely surjective.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -38, 90])  # 361A
            sage: E.galois_representation().non_surjective()   # CM curve
            [0]

        ::

            sage: E = EllipticCurve([0, -1, 1, 0, 0]) # X_1(11)
            sage: E.galois_representation().non_surjective()
            [5]

        ::

            sage: E = EllipticCurve([0, 0, 1, -1, 0]) # 37A
            sage: E.galois_representation().non_surjective()
            []

        ::

            sage: E = EllipticCurve([0,-1,1,-2,-1])   # 141C
            sage: E.galois_representation().non_surjective()
            [13]

        ALGORITHM: When `p = 3` use division polynomials.
        For `5 = p = B`, where
        `B` is Cojocaru's bound, use the results in Section 2 of

        """
        if self.E.has_cm():
            misc.verbose("cm curve")
            return [0]
        N = self.E.conductor()
        if self.E.is_semistable():
            C = 11
            misc.verbose("semistable -- so bound is 11")
        else:
            C = 1 + 4*sqrt(6)*int(N)/3 * sqrt(misc.mul([1+1.0/int(p) for p,_ in arith.factor(N)]))
            misc.verbose("conductor = %s, and bound is %s"%(N,C))
        C = 1 + 4*sqrt(6)*int(N)/3 * sqrt(misc.mul([1+1.0/int(p) for p,_ in arith.factor(N)]))
        misc.verbose("conductor = %s, and bound is %s"%(N,C))
        B = []
        p = 2
        while p <= C:
            t = self.is_surjective(p, A=A)
            misc.verbose("(%s,%s)"%(p,t))
            if not t:
                B.append(p)
            p = arith.next_prime(p)
        return B
