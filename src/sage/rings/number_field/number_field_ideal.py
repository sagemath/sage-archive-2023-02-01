"""
Number Field Ideals

AUTHOR:
   -- Steven Sivek (2005-05-16)
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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

import operator

import sage.rings.field_element as field_element
import sage.rings.polynomial_element as polynomial
import sage.rings.polynomial_ring as polynomial_ring
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.arith as arith
import sage.misc.misc as misc

import number_field
import number_field_element

from sage.libs.all import pari_gen
from sage.rings.ideal import (Ideal_generic, Ideal_fractional)
from sage.misc.misc import prod

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()

def is_NumberFieldIdeal(x):
    return isinstance(x, NumberFieldIdeal)

def convert_from_zk_basis(field, hnf):
    return field.pari_nf().getattr('zk') * hnf

class NumberFieldIdeal(Ideal_fractional):
    """
    An ideal of a number field.
    """
    def __init__(self, field, gens, coerce=True):
        """
        INPUT:
            field -- a number field
            x -- a list of NumberFieldElements belonging to the field
        """
        if not isinstance(field, number_field.NumberField_generic):
            raise TypeError, "field (=%s) must be a number field."%field

        Ideal_generic.__init__(self, field, gens, coerce)

    def __cmp__(self, other):
        if self.pari_hnf() == other.pari_hnf():
            return 0
        return -1

    def _coerce_impl(self, x):
        if isinstance(x, NumberFieldIdeal):
            if x.parent() == self.parent():
                return x
        if isinstance(x, (rational.Rational, integer.Integer, int, long)):
            return self.number_field().ideal(x)
        if isinstance(x, NumberFieldElement) and x.parent() == self.number_field():
            return self.number_field().ideal(x)
        if isinstance(x, (list, tuple)):
            return self.number_field().ideal(list(x))
        raise TypeError

    def _contains_(self, x):
        """
        For now, $x \in I$ if and only if $\langle x \rangle + I = I$.
        """
        x_ideal = self.number_field().ideal(x)
        return self + x_ideal == self

    def __elements_from_hnf(self, hnf):
        """
        Convert a matrix in Hermite normal form (from PARI) to a list of
        NumberFieldElements.
        """
        K = self.number_field()
        nf = K.pari_nf()
        R = K.polynomial().parent()
        return [ K(R(x)) for x in convert_from_zk_basis(K, hnf) ]

    def _repr_short(self):
        return '(%s)'%(', '.join([str(x) for x in self.gens_reduced()]))

    def _repr_(self):
        if self.is_integral():
            return Ideal_generic._repr_(self)
        else:
            return self.__repr__()

    def __div__(self, other):
        """
        Return the quotient self / other.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2 - 5)
            sage: I = K.ideal(2/(5+a))
            sage: J = K.ideal(17+a)
            sage: I/J
            Fractional ideal (-17/1420*a + 1/284) of Number Field in a with defining polynomial x^2 - 5
            sage: (I/J) * J
            Fractional ideal (-1/5*a) of Number Field in a with defining polynomial x^2 - 5
            sage: (I/J) * J == I
            True
        """
        return self * other.__invert__()

    def __invert__(self):
        """
        Return the multiplicative inverse of self.  Call with ~self.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: I = K.ideal(2/(5+a))
            sage: ~I
            Fractional ideal (1/2*a + 5/2) of Number Field in a with defining polynomial x^3 - 2
            sage: 1/I
            Fractional ideal (1/2*a + 5/2) of Number Field in a with defining polynomial x^3 - 2
            sage: (1/I) * I
            Fractional ideal (1) of Number Field in a with defining polynomial x^3 - 2
        """
        if self.is_zero():
            raise ZeroDivisionError
        nf = self.number_field().pari_nf()
        hnf = nf.idealdiv(self.number_field().ideal(1).pari_hnf(),
                          self.pari_hnf())
        I = self.number_field().ideal(self.__elements_from_hnf(hnf))
        I.__pari_hnf = hnf
        return I

    def __pow__(self, right):
        """
        Return self to the power of right.
        """
        right = int(right)
        if right < 0:
            x = self.number_field().ideal(1) / self
            right *= -1
            return arith.generic_power(x, right, one=self.parent()(1))
        return arith.generic_power(self, right, one=self.parent()(1))

    def pari_hnf(self):
        """
        Return PARI's representation of this ideal in Hermite normal form.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: I = K.ideal(2/(5+a))
            sage: I.pari_hnf()
            [2, 0, 50/127; 0, 2, 244/127; 0, 0, 2/127]
        """
        try:
            return self.__pari_hnf
        except AttributeError:
            nf = self.number_field().pari_nf()
            self.__pari_hnf = nf.idealhnf(0)
            hnflist = [ nf.idealhnf(x.polynomial()) for x in self.gens() ]
            for ideal in hnflist:
                self.__pari_hnf = nf.idealadd(self.__pari_hnf, ideal)
            return self.__pari_hnf

    def divides(self, other):
        """
        Returns True if this ideal divides other and False otherwise.
        """
        if not isinstance(other, NumberFieldIdeal):
            other = self.number_field().ideal(other)
        if self.is_zero():
            return other.is_zero # since 0 \subset 0
        return (other / self).is_integral()

    def factor(self):
        """
        Factorization of this ideal in terms of prime ideals.
        Output is of the form [ (ideal, exponent), (ideal, exponent), ... ].
        """
        try:
            return self.__factorization
        except AttributeError:
            if self.is_zero():
                self.__factorization = []
                return self.__factorization
            K = self.number_field()
            F = list(K.pari_nf().idealfactor(self.pari_hnf()))
            P, exps = F[0], F[1]
            A = []
            zk_basis = K.pari_nf().getattr('zk')
            for i, p in enumerate(P):
                prime, alpha = p.getattr('gen')
                I = K.ideal([ZZ(prime), K(zk_basis * alpha)])
                I.__pari_prime = p
                A.append((I,ZZ(exps[i])))
            self.__factorization = A
            return self.__factorization

    def gens_reduced(self, certify=True):
        """
        Express this ideal in terms of at most two generators, and
        one if possible.  Note that if the ideal is not principal,
        then this uses PARI's \code{idealtwoelt} function, which
        takes exponential time, the first time it is called for each
        ideal.  Also, this indirectly uses \code{bnfisprincipal}, so
        set \code{certify=False} if you don't want to prove correctness.

        EXAMPLE:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1, 'i')
            sage: J = K.ideal([i+1, 2])
            sage: J.gens()
            (2, i + 1)
            sage: J.gens_reduced()
            (i + 1,)
        """
        try:
            ## Compute the single generator, if it exists
            dummy = self.is_principal(certify)
            return self.__reduced_generators
        except AttributeError:
            K = self.number_field()
            nf = K.pari_nf()
            R = K.polynomial().parent()
            if self.is_prime():
                a = self.smallest_integer()
                alpha = nf.idealtwoelt(self.pari_hnf(), a)
            else:
                a, alpha = nf.idealtwoelt(self.pari_hnf())
            gens = [ QQ(a), K(R(nf.getattr('zk')*alpha)) ]
            if gens[1] in self.number_field().ideal(gens[0]):
                gens = [ gens[0] ]
            elif gens[0] in self.number_field().ideal(gens[1]):
                gens = [ gens[1] ]
            self.__reduced_generators = tuple(gens)
            return self.__reduced_generators

    def integral_basis(self):
        r"""
        Return a list of generators for this ideal as a $\mathbb{Z}$-module.

        EXAMPLE:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)
            sage: J = K.ideal(i+1)
            sage: J.integral_basis()
            [2, i + 1]
        """
        hnf = self.pari_hnf()
        return self.__elements_from_hnf(hnf)

    def integral_split(self):
        """
        Return a tuple (I, d), where I is an integral ideal, and d is the
        smallest positive integer such that this ideal is equal to I/d.

        EXAMPLE:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2-5)
            sage: I = K.ideal(2/(5+a))
            sage: I.is_integral()
            False
            sage: J,d = I.integral_split()
            sage: J
            Fractional ideal (-1/2*a + 5/2) of Number Field in a with defining polynomial x^2 - 5
            sage: J.is_integral()
            True
            sage: d
            5
            sage: I == J/d
            True
        """
        try:
            return self.__integral_split
        except AttributeError:
            if self.is_integral():
                self.__integral_split = (self, ZZ(1))
            else:
                factors = self.factor()
                denom_list = filter(lambda (p,e): e < 0 , factors)
                denominator = prod([ p.smallest_integer()**(-e)
                                     for (p,e) in denom_list ])
                ## Get a list of the primes dividing the denominator
                plist = [ p.smallest_integer() for (p,e) in denom_list ]
                for p in plist:
                    while denominator % p == 0 and (self*(denominator/p)).is_integral():
                        denominator //= p
                self.__integral_split = (self*denominator, denominator)
            return self.__integral_split

    def is_integral(self):
        """
        Determine if this ideal is integral.

        EXAMPLES:
           sage: R.<x> = PolynomialRing(QQ)
           sage: K.<a> = NumberField(x^5-x+1)
           sage: K.ideal(a).is_integral()
           True
           sage: (K.ideal(1) / (3*a+1)).is_integral()
           False
        """
        try:
            return self.__is_integral
        except AttributeError:
            one = self.number_field().ideal(1)
            self.__is_integral = (self+one == one)
            return self.__is_integral

    def is_maximal(self):
        """
        Determine if this ideal is maximal.
        """
        return self.is_prime() and not self.is_zero()

    def is_prime(self):
        """
        Determine if this ideal is prime.
        """
        try:
            return self.__pari_prime is not None
        except AttributeError:
            if self.is_zero():
                self.__pari_prime = []
                return True
            K = self.number_field()
            F = list(K.pari_nf().idealfactor(self.pari_hnf()))
            P, exps = F[0], F[1]
            if len(P) != 1 or exps[0] != 1:
                self.__pari_prime = None
            else:
                self.__pari_prime = P[0]
            return self.__pari_prime is not None

    def is_principal(self, certify=True):
        """
        Determine whether this ideal is principal.  Since it uses the
        PARI method \code{bnfisprincipal}, specify \code{certify=True}
        (the default setting) to prove the correctness of the output.
        """
        try:
            return self.__is_principal
        except AttributeError:
            if len (self.gens()) <= 1:
                self.__is_principal = True
                self.__reduced_generators = self.gens()
                return self.__is_principal
            bnf = self.number_field().pari_bnf(certify)
            v = bnf.bnfisprincipal(self.pari_hnf())
            self.__is_principal = (v[0] == 0) ## i.e., v[0] is the zero vector
            if self.__is_principal:
                K = self.number_field()
                R = K.polynomial().parent()
                g = K(R(bnf.getattr('zk') * v[1]))
                self.__reduced_generators = tuple([g])
            return self.__is_principal

    def is_zero(self):
        """
        Determine if this is the zero ideal.
        """
        return self == self.number_field().ideal(0)

    def norm(self):
        return QQ(self.number_field().pari_nf().idealnorm(self.pari_hnf()))

    def number_field(self):
        return self.ring()

    def ramification(self):
        r"""
        Return the ramification index of this ideal, assuming it is prime
        and nonzero.  Otherwise, raise a ValueError.
        """
        if self.is_zero():
            raise ValueError, "The ideal (=%s) is zero"%self
        if self.is_prime():
            return ZZ(self.__pari_prime.getattr('e'))
        raise ValueError, "the ideal (= %s) is not prime"%self

    def residue_class_degree(self):
        r"""
        Return the residue class degree of this ideal, assuming it is prime
        and nonzero.  Otherwise, raise a ValueError.
        """
        if self.is_zero():
            raise ValueError, "The ideal (=%s) is zero"%self
        if self.is_prime():
            return ZZ(self.__pari_prime.getattr('f'))
        raise ValueError, "the ideal (= %s) is not prime"%self

    def smallest_integer(self):
        r"""
        Return the smallest nonnegative integer in $I \cap \mathbb{Z}$,
        where $I$ is this ideal.  If $I = 0$, raise a ValueError.

        EXAMPLE:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2+6)
            sage: I = K.ideal([4,a])/7
            sage: I.smallest_integer()
            2
        """
        if self.is_zero():
            raise ValueError, "ideal (= %s) must be nonzero"%self
        try:
            return self.__smallest_integer
        except AttributeError:
            if self.is_prime():
                self.__smallest_integer = ZZ(self.__pari_prime.getattr('p'))
                return self.__smallest_integer
            if self.is_integral():
                factors = self.factor()
                bound = prod([ p.smallest_integer()**e for (p,e) in factors ])
                plist = [ p.smallest_integer() for (p,e) in factors ]
                plist.sort()
                indices = filter(lambda(i): i==0 or plist[i] != plist[i-1],
                                 range(0,len(plist)))
                plist = [ plist[i] for i in indices ] ## unique list of primes
                for p in plist:
                    while bound % p == 0 and (self/(bound/p)).is_integral():
                        bound /= p
                self.smallest_integer = ZZ(bound)
                return self.__smallest_integer
            I,d = self.integral_split() ## self = I/d
            n = I.smallest_integer()    ## n/d in self
            self.__smallest_integer =  n / arith.gcd(ZZ(n),ZZ(d))
            return self.__smallest_integer

    def valuation(self, p):
        r"""
        Return the valuation of this ideal at the prime $\mathfrak{p}$.
        If $\mathfrak{p}$ is not prime, raise a ValueError.
        """
        if not isinstance(p, NumberFieldIdeal):
            p = self.number_field().ideal(p)
        if p.is_zero() or not p.is_prime():
            raise ValueError, "p (= %s) must be a nonzero prime"%p
        nf = self.number_field().pari_nf()
        return ZZ(nf.idealval(self.pari_hnf(), p.__pari_prime))



class NumberFieldIdeal_rel(NumberFieldIdeal):
    """
    An ideal of a relative number field.
    """
    def pari_rhnf(self):
        """
        Return PARI's representation of this relative ideal in Hermite
        normal form.
        """
        try:
            return self.__pari_rhnf
        except AttributeError:
            nf = self.number_field().absolute_field().pari_nf()
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
            L = self.number_field().absolute_field()
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
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2+6)
            sage: L.<b> = K.extension(K['x'].gen()^4 + a)
            sage: L.ideal(b).norm()
            Fractional ideal (-a) of Number Field in a with defining polynomial x^2 + 6
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
