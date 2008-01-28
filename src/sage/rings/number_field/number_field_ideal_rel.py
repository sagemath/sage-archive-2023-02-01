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

class NumberFieldIdeal_rel(NumberFieldFractionalIdeal):
    """
    An ideal of a relative number field.

    EXAMPLES:
        sage: K.<a> = NumberField([x^2 + 1, x^2 + 2]); K
        Number Field in a0 with defining polynomial x^2 + 1 over its base field
        sage: i = K.ideal(38); i
        Fractional ideal (38)

        sage: K.<a> = NumberField([x^2 + 1, x^2 + 2]); K
        Number Field in a0 with defining polynomial x^2 + 1 over its base field
        sage: i = K.ideal([a+1]); i
        Fractional ideal (a0 + 1)
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
            rnf = self.number_field().pari_rnf()
            L = self.number_field().absolute_field('a')
            R = L['x']
            nf = L.pari_nf()
            genlist = [L(R(x.polynomial())) for x in list(self.gens())]
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
