"""
Relative Number Field Ideals

AUTHOR:
   -- Steven Sivek (2005-05-16)
   -- Willia Stein (2007-09-06)
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
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
#*****************************************************************************

from number_field_ideal import NumberFieldFractionalIdeal, convert_from_zk_basis

import sage.rings.rational_field as rational_field
QQ = rational_field.RationalField()

class NumberFieldFractionalIdeal_rel(NumberFieldFractionalIdeal):
    """
    An ideal of a relative number field.

    EXAMPLES:
        sage: K.<a> = NumberField([x^2 + 1, x^2 + 2]); K
        Number Field in a0 with defining polynomial x^2 + 1 over its base field
        sage: i = K.ideal(38); i
        Fractional ideal (38)

    WARNING: Ideals in relative number fields are broken:
        sage: K.<a> = NumberField([x^2 + 1, x^2 + 2]); K
        Number Field in a0 with defining polynomial x^2 + 1 over its base field
        sage: i = K.ideal([a+1]); i
        Traceback (most recent call last):
        ...
        TypeError: Unable to coerce -a1 to a rational
    """
    def pari_rhnf(self):
        """
        Return PARI's representation of this relative ideal in Hermite
        normal form.

        EXAMPLES:

        """
        try:
            return self.__pari_rhnf
        except AttributeError:
            nf = self.number_field().absolute_field('a').pari_nf()
            rnf = self.number_field().pari_rnf()
            L_hnf = self.absolute_ideal().pari_hnf()
            self.__pari_rhnf = rnf.rnfidealabstorel(nf.getattr('zk')*L_hnf)
            return self.__pari_rhnf

    def absolute_ideal(self):
        """
        If this is an ideal in the extension L/K, return the ideal with
        the same generators in the absolute field L/Q.
        """
        try:
            return self.__absolute_ideal
        except AttributeError:
            L = self.number_field().absolute_field('a')
            genlist = [L(x.polynomial()) for x in list(self.gens())]
            self.__absolute_ideal = L.ideal(genlist)
            return self.__absolute_ideal

    def gens_reduced(self):
        try:
            return self.__reduced_generators
        except AttributeError:
            L = self.number_field()
            K = L.base_field()
            R = L.polynomial().parent()
            S = L['x']
            gens = L.pari_rnf().rnfidealtwoelt(self.pari_rhnf())
            gens = [ L(R(x.lift().lift())) for x in gens ]
            ## Make sure that gens[1] is in L, not K
            Lcoeff = [ L(x) for x in list(gens[1].polynomial()) ]
            gens[1] = S.hom([L.gen()])(S(Lcoeff))

            if gens[1] in L.ideal([ gens[0] ]):
                gens = [ gens[0] ]
            elif gens[0] in L.ideal([ gens[1] ]):
                gens = [ gens[1] ]
            self.__reduced_generators = tuple(gens)
            return self.__reduced_generators

    def __invert__(self):
        """
        Return the multiplicative inverse of self.  Call with ~self.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 + 2])
            sage: I = K.fractional_ideal(4)
            sage: I^(-1)
            Fractional ideal (1/4)
            sage: I * I^(-1)
            Fractional ideal (1)
        """
        if self.is_zero():
            raise ZeroDivisionError
        L = self.number_field()
        R = QQ['x']
        rnf = L.pari_rnf()
        inverse = self.absolute_ideal().__invert__()
        nf_zk = inverse.number_field().pari_nf().getattr('zk')
        genlist = [L(R(x)) for x in nf_zk * inverse.pari_hnf()]
        return L.ideal(genlist)

    def is_principal(self):
        return self.absolute_ideal().is_principal()

    def is_zero(self):
        zero = self.number_field().pari_rnf().rnfidealhnf(0)
        return self.pari_rhnf() == zero

    def norm(self):
        """
        Compute the relative norm of this extension L/K as an ideal of K.

        EXAMPLE:
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2+6)
            sage: L.<b> = K.extension(K['x'].gen()^4 + a)
            sage: N = L.ideal(b).norm(); N
            Fractional ideal (-a)
            sage: N.parent()
            Monoid of ideals of Number Field in a with defining polynomial x^2 + 6
            sage: N.ring()
            Number Field in a with defining polynomial x^2 + 6
        """
        L = self.number_field()
        K = L.base_field()
        R = K.polynomial().parent()
        hnf = L.pari_rnf().rnfidealnormrel(self.pari_rhnf())
        return K.ideal([ K(R(x)) for x in convert_from_zk_basis(K, hnf) ])

    def ideal_below(self):
        """
        Compute the ideal of K below this ideal of L.

        EXAMPLE:
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2+6)
            sage: L.<b> = K.extension(K['x'].gen()^4 + a)
            sage: N = L.ideal(b)
            sage: M = N.ideal_below(); M == K.ideal([-a])
            True
            sage: Np = L.ideal( [ L(t) for t in M.gens() ])
            sage: Np.ideal_below() == M
            True
            sage: M.parent()
            Monoid of ideals of Number Field in a with defining polynomial x^2 + 6
            sage: M.ring()
            Number Field in a with defining polynomial x^2 + 6
            sage: M.ring() is K
            True

            This example concerns an inert ideal:

            sage: K = NumberField(x^4 + 6*x^2 + 24, 'a')
            sage: K.factor(7)
            Fractional ideal (7)
            sage: K0, K0_into_K, _ = K.subfields(2)[0]
            sage: K0
            Number Field in a0 with defining polynomial x^2 - 6*x + 24
            sage: L = K.relativize(K0_into_K, 'c'); L
            Number Field in c0 with defining polynomial x^2 + a0 over its base field
            sage: L.base_field() is K0
            True
            sage: L.ideal(7)
            Fractional ideal (7)
            sage: L.ideal(7).ideal_below()
            Fractional ideal (7)
            sage: L.ideal(7).ideal_below().number_field() is K0
            True

            This example concerns an ideal that splits in the quadratic field
            but each factor ideal remains inert in the extension:

            sage: len(K.factor(19))
            2
            sage: K0 = L.base_field(); a0 = K0.gen()
            sage: len(K0.factor(19))
            2
            sage: w1 = -a0 + 1; P1 = K0.ideal([w1])
            sage: P1.norm().factor(), P1.is_prime()
            (19, True)
            sage: L_into_K, K_into_L = L.structure()
            sage: L.ideal(K_into_L(K0_into_K(w1))).ideal_below() == P1
            True

            The choice of embedding of quadratic field into quartic field
            matters:

            sage: rho, tau = K0.embeddings(K)
            sage: L1 = K.relativize(rho, 'b')
            sage: L2 = K.relativize(tau, 'b')
            sage: L1_into_K, K_into_L1 = L1.structure()
            sage: L2_into_K, K_into_L2 = L2.structure()
            sage: a = K.gen()
            sage: P = K.ideal([a^2 + 5])
            sage: K_into_L1(P).ideal_below() == K0.ideal([-a0 + 1])
            True
            sage: K_into_L2(P).ideal_below() == K0.ideal([-a0 + 5])
            True
            sage: K0.ideal([-a0 + 1]) == K0.ideal([-a0 + 5])
            False
        """
        L = self.number_field()
        K = L.base_field()
        R = K.polynomial().parent()
        hnf = L.pari_rnf().rnfidealdown(self.pari_rhnf())
        return K.ideal([ K(R(x)) for x in convert_from_zk_basis(K, hnf) ])

    def factor(self):
        raise NotImplementedError
    def integral_basis(self):
        raise NotImplementedError
    def integral_split(self):
        raise NotImplementedError
    def is_maximal(self):
        raise NotImplementedError
    def is_prime(self):
        raise NotImplementedError
    def ramification(self):
        raise NotImplementedError
    def residue_class_degree(self):
        raise NotImplementedError
    def smallest_integer(self):
        raise NotImplementedError
    def valuation(self):
        raise NotImplementedError

def is_NumberFieldFractionalIdeal_rel(x):
    """
    Return True if x is a fractional ideal of a relative number field.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field_ideal_rel import is_NumberFieldFractionalIdeal_rel
        sage: from sage.rings.number_field.number_field_ideal import is_NumberFieldFractionalIdeal
        sage: is_NumberFieldFractionalIdeal_rel(2/3)
        False
        sage: is_NumberFieldFractionalIdeal_rel(ideal(5))
        False
        sage: k.<a> = NumberField(x^2 + 2)
        sage: I = k.ideal([a + 1]); I
        Fractional ideal (a + 1)
        sage: is_NumberFieldFractionalIdeal_rel(I)
        False
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^2+6)
        sage: L.<b> = K.extension(K['x'].gen()^4 + a)
        sage: I = L.ideal(b); I
        Fractional ideal (b)
        sage: is_NumberFieldFractionalIdeal_rel(I)
        True
        sage: N = I.norm(); N
        Fractional ideal (-a)
        sage: is_NumberFieldFractionalIdeal_rel(N)
        False
        sage: is_NumberFieldFractionalIdeal(N)
        True
    """
    return isinstance(x, NumberFieldFractionalIdeal_rel)

