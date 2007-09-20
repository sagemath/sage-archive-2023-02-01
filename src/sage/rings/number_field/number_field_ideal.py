"""
Number Field Ideals

AUTHOR:
   -- Steven Sivek (2005-05-16)
   -- Willia Stein (2007-09-06): vastly improved the doctesting

TESTS:
We test that pickling works:
    sage: K.<a> = NumberField(x^2 - 5)
    sage: I = K.ideal(2/(5+a))
    sage: I == loads(dumps(I))
    True
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
import sage.rings.polynomial.polynomial_element as polynomial
import sage.rings.polynomial.polynomial_ring as polynomial_ring
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.rational as rational
import sage.rings.integer as integer
import sage.rings.arith as arith
import sage.misc.misc as misc

import number_field
import number_field_element

from sage.libs.all import pari_gen
from sage.rings.ideal import (Ideal_generic, Ideal_fractional)
from sage.misc.misc import prod
from sage.structure.element import generic_power
from sage.structure.factorization import Factorization

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()

def is_NumberFieldIdeal(x):
    """
    Return True if x is a fractional ideal of a number field.

    EXAMPLES:
        sage: is_NumberFieldIdeal(2/3)
        False
        sage: is_NumberFieldIdeal(ideal(5))
        False
        sage: k.<a> = NumberField(x^2 + 2)
        sage: I = k.ideal([a + 1]); I
        Fractional ideal (a + 1) of Number Field in a with defining polynomial x^2 + 2
        sage: is_NumberFieldIdeal(I)
        True
    """
    return isinstance(x, NumberFieldIdeal)

def convert_from_zk_basis(field, hnf):
    """
    Used internally in the number field ideal implementation for
    converting from the order basis to the number field basis.

    INPUT:
        field -- a number field
        hnf -- a pari HNF matrix, output by the pari_hnf() function.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field_ideal import convert_from_zk_basis
        sage: k.<a> = NumberField(x^2 + 23)
        sage: I = k.factor_integer(3)[0][0]; I
        Fractional ideal (3, -1/2*a + 1/2) of Number Field in a with defining polynomial x^2 + 23
        sage: hnf = I.pari_hnf(); hnf
        [3, 0; 0, 1]
        sage: convert_from_zk_basis(k, hnf)
        [3, 1/2*x - 1/2]

    """
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

        EXAMPLES:
            sage: NumberField(x^2 + 1, 'a').ideal(7)
            Fractional ideal (7) of Number Field in a with defining polynomial x^2 + 1
        """
        if not isinstance(field, number_field.NumberField_generic):
            raise TypeError, "field (=%s) must be a number field."%field

        Ideal_generic.__init__(self, field, gens, coerce)

    def __cmp__(self, other):
        """
        Compare these a fractional ideal of a number field to
        something else.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 3); K
            Number Field in a with defining polynomial x^2 + 3
            sage: f = K.factor_integer(15); f
            (Fractional ideal (1/2*a - 3/2) of Number Field in a with defining polynomial x^2 + 3)^2 * (Fractional ideal (5) of Number Field in a with defining polynomial x^2 + 3)
            sage: cmp(f[0][0], f[1][0])
            -1
            sage: cmp(f[0][0], f[0][0])
            0
            sage: cmp(f[1][0], f[0][0])
            1
            sage: f[1][0] == 5
            True
            sage: f[1][0] == GF(7)(5)
            False
        """
        if not isinstance(other, NumberFieldIdeal):
            return cmp(type(self), type(other))
        return cmp(self.pari_hnf(), other.pari_hnf())

    def _contains_(self, x):
        """
        Return True if x is an element of this fractional ideal.

        This function is called (indirectly) when the \code{in}
        operator is used.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 23); K
            Number Field in a with defining polynomial x^2 + 23
            sage: I = K.factor_integer(13)[0][0]; I
            Fractional ideal (13, a - 4) of Number Field in a with defining polynomial x^2 + 23
            sage: a in I
            False
            sage: 13 in I
            True
            sage: 13/2 in I
            False
            sage: a + 9 in I
            True

            sage: K.<a> = NumberField(x^4 + 3); K
            Number Field in a with defining polynomial x^4 + 3
            sage: I = K.factor_integer(13)[0][0]
            sage: I  # random sign in output
            Fractional ideal (-2*a^2 - 1) of Number Field in a with defining polynomial x^4 + 3
            sage: 2/3 in I
            False
            sage: 1 in I
            False
            sage: 13 in I
            True
            sage: 1 in I*I^(-1)
            True
            sage: I   # random sign in output
            Fractional ideal (-2*a^2 - 1) of Number Field in a with defining polynomial x^4 + 3
        """
        # For now, $x \in I$ if and only if $\langle x \rangle + I = I$.
        # Is there a better way to do this?
        x_ideal = self.number_field().ideal(x)
        return self + x_ideal == self

    def __elements_from_hnf(self, hnf):
        """
        Convert a PARI Hermite normal form matrix to a list of
        NumberFieldElements.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 389); K
            Number Field in a with defining polynomial x^3 + 389
            sage: I = K.factor_integer(17)[0][0]; I
            Fractional ideal (-100*a^2 + 730*a - 5329) of Number Field in a with defining polynomial x^3 + 389
            sage: hnf = I.pari_hnf(); hnf
            [17, 0, 13; 0, 17, 8; 0, 0, 1]
            sage: I._NumberFieldIdeal__elements_from_hnf(hnf)
            [17, 17*a, a^2 + 8*a + 13]
            sage: I._NumberFieldIdeal__elements_from_hnf(hnf^(-1))
            [1/17, 1/17*a, a^2 - 8/17*a - 13/17]
        """
        K = self.number_field()
        nf = K.pari_nf()
        R = K.polynomial().parent()
        return [ K(R(x)) for x in convert_from_zk_basis(K, hnf) ]

    def _repr_short(self):
        """
        Efficient string representation of this fraction ideal.

        EXAMPLES:
            sage: K.<a> = NumberField(x^4 + 389); K
            Number Field in a with defining polynomial x^4 + 389
            sage: I = K.factor_integer(17)[0][0]; I
            Fractional ideal (17, a^2 - 6) of Number Field in a with defining polynomial x^4 + 389
            sage: I._repr_short()
            '(17, a^2 - 6)'
        """
        return '(%s)'%(', '.join([str(x) for x in self.gens_reduced()]))

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

    def __pow__(self, r):
        """
        Return self to the power of right.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: I = K.ideal(2/(5+a))
            sage: J = I^2
            sage: K = I^(-2)
            sage: J*K
            Fractional ideal (1) of Number Field in a with defining polynomial x^3 - 2
        """
        return generic_power(self, r)

    def _pari_(self):
        return self.pari_hnf()

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

        EXAMPLES:
            sage: K.<a> = CyclotomicField(11); K
            Cyclotomic Field of order 11 and degree 10
            sage: I = K.factor_integer(31)[0][0]; I
            Fractional ideal (-3*a^7 - 4*a^5 - 3*a^4 - 3*a^2 - 3*a - 3) of Cyclotomic Field of order 11 and degree 10
            sage: I.divides(I)
            True
            sage: I.divides(31)
            True
            sage: I.divides(29)
            False
        """
        if not isinstance(other, NumberFieldIdeal):
            other = self.number_field().ideal(other)
        if self.is_zero():
            return other.is_zero # since 0 \subset 0
        return (other / self).is_integral()

    def factor(self):
        """
        Factorization of this ideal in terms of prime ideals.

        EXAMPLES:
            sage: K.<a> = NumberField(x^4 + 23); K
            Number Field in a with defining polynomial x^4 + 23
            sage: I = K.ideal(19); I
            Fractional ideal (19) of Number Field in a with defining polynomial x^4 + 23
            sage: F = I.factor(); F
            (Fractional ideal (a^2 + 2*a + 2) of Number Field in a with defining polynomial x^4 + 23) * (Fractional ideal (a^2 - 2*a + 2) of Number Field in a with defining polynomial x^4 + 23)
            sage: type(F)
            <class 'sage.structure.factorization.Factorization'>
            sage: list(F)
            [(Fractional ideal (a^2 + 2*a + 2) of Number Field in a with defining polynomial x^4 + 23, 1), (Fractional ideal (a^2 - 2*a + 2) of Number Field in a with defining polynomial x^4 + 23, 1)]
            sage: F.prod()
            Fractional ideal (19) of Number Field in a with defining polynomial x^4 + 23
        """
        try:
            return self.__factorization
        except AttributeError:
            if self.is_zero():
                self.__factorization = Factorization([])
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
            self.__factorization = Factorization(A)
            return self.__factorization

    def reduce_equiv(self):
        """
        Return a small ideal that is equivalent to self in the group
        of fractional ideals modulo principal ideals.  Very often (but
        not always) if self is principal then this function returns
        the unit ideal.

        ALGORITHM: Calls pari's idealred function.
        """
        K = self.number_field()
        P = K.pari_nf()
        hnf = P.idealred(self.pari_hnf())
        gens = self.__elements_from_hnf(hnf)
        return K.ideal(gens)

    def gens_reduced(self, proof=None):
        r"""
        Express this ideal in terms of at most two generators, and one
        if possible.

        Note that if the ideal is not principal, then this uses PARI's
        \code{idealtwoelt} function, which takes exponential time, the
        first time it is called for each ideal.  Also, this indirectly
        uses \code{bnfisprincipal}, so set \code{proof=True} if you
        want to prove correctness (which \emph{is} the default).

        EXAMPLE:
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1, 'i')
            sage: J = K.ideal([i+1, 2])
            sage: J.gens()
            (i + 1, 2)
            sage: J.gens_reduced()
            (i + 1,)
        """
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "number_field")
        try:
            ## Compute the single generator, if it exists
            dummy = self.is_principal(proof)
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
        Return True if this ideal is integral.

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
        Return True if this ideal is maximal.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 3); K
            Number Field in a with defining polynomial x^3 + 3
            sage: K.ideal(5).is_maximal()
            False
            sage: K.ideal(7).is_maximal()
            True
        """
        return self.is_prime() and not self.is_zero()

    def is_prime(self):
        """
        Return True if this ideal is prime.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 - 17); K
            Number Field in a with defining polynomial x^2 - 17
            sage: K.ideal(5).is_prime()
            True
            sage: K.ideal(13).is_prime()
            False
            sage: K.ideal(17).is_prime()
            False
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

    def is_principal(self, proof=None):
        r"""
        Return True if this ideal is principal.

        Since it uses the PARI method \code{bnfisprincipal}, specify
        \code{proof=True} (this is the default setting) to prove the
        correctness of the output.

        EXAMPLES:
        We create equal ideals in two different ways, and note that
        they are both actually principal ideals.
            sage: K = QuadraticField(-119,'a')
            sage: P = K.ideal([2]).factor()[1][0]
            sage: I = P^5
            sage: I.is_principal()
            True
        """
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "number_field")
        try:
            return self.__is_principal
        except AttributeError:
            if len (self.gens()) <= 1:
                self.__is_principal = True
                self.__reduced_generators = self.gens()
                return self.__is_principal
            bnf = self.number_field().pari_bnf(proof)
            v = bnf.bnfisprincipal(self.pari_hnf())
            self.__is_principal = is_pari_zero_vector(v[0])
            if self.__is_principal:
                K = self.number_field()
                R = K.polynomial().parent()
                g = K(R(bnf.getattr('zk') * v[1]))
                self.__reduced_generators = tuple([g])
            return self.__is_principal

    def is_zero(self):
        """
        Return True if this is the zero ideal.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 2); K
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(3).is_zero()
            False
            sage: K.ideal(0).is_zero()
            True
        """
        return self == self.number_field().ideal(0)

    def norm(self):
        """
        Return the norm of this fractional ideal as a rational number.

        EXAMPLES:
            sage: K.<a> = NumberField(x^4 + 23); K
            Number Field in a with defining polynomial x^4 + 23
            sage: I = K.ideal(19); I
            Fractional ideal (19) of Number Field in a with defining polynomial x^4 + 23
            sage: factor(I.norm())
            19^4
            sage: F = I.factor()
            sage: F[0][0].norm().factor()
            19^2
        """
        return QQ(self.number_field().pari_nf().idealnorm(self.pari_hnf()))

    def number_field(self):
        """
        Return the number field that this is a fractional ideal in.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 2); K
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(3).number_field()
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(0).number_field()
            Number Field in a with defining polynomial x^2 + 2
        """
        return self.ring()

    def ramification_index(self):
        r"""
        Return the ramification index of this ideal, assuming it is prime
        and nonzero.  Otherwise, raise a ValueError.

        The ramification index is the power of this prime appearing in
        the factorization of the prime in $\ZZ$ that this primes lies
        over.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 2); K
            Number Field in a with defining polynomial x^2 + 2
            sage: f = K.factor_integer(2); f
            (Fractional ideal (-a) of Number Field in a with defining polynomial x^2 + 2)^2
            sage: f[0][0].ramification_index()
            2
            sage: K.ideal(13).ramification_index()
            1
            sage: K.ideal(17).ramification_index()
            Traceback (most recent call last):
            ...
            ValueError: the ideal (= Fractional ideal (17) of Number Field in a with defining polynomial x^2 + 2) is not prime
        """
        if self.is_zero():
            raise ValueError, "The ideal (=%s) is zero"%self
        if self.is_prime():
            return ZZ(self.__pari_prime.getattr('e'))
        raise ValueError, "the ideal (= %s) is not prime"%self

    def residue_class_degree(self):
        r"""
        Return the residue class degree of this ideal, assuming it is
        prime and nonzero.  Otherwise, raise a ValueError.

        The residue class degree of a prime ideal $I$ is the degree of
        the extension $O_K/I$ of its prime subfield.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 + 2); K
            Number Field in a with defining polynomial x^5 + 2
            sage: f = K.factor_integer(19); f
            (Fractional ideal (a^2 + a - 3) of Number Field in a with defining polynomial x^5 + 2) * (Fractional ideal (-2*a^4 - a^2 + 2*a - 1) of Number Field in a with defining polynomial x^5 + 2) * (Fractional ideal (a^2 + a - 1) of Number Field in a with defining polynomial x^5 + 2)
            sage: [i.residue_class_degree() for i, _ in f]
            [2, 2, 1]
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

        INPUT:
            p -- a prime ideal of this number field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 + 2); K
            Number Field in a with defining polynomial x^5 + 2
            sage: i = K.ideal(38); i
            Fractional ideal (38) of Number Field in a with defining polynomial x^5 + 2
            sage: i.valuation(K.factor_integer(19)[0][0])
            1
            sage: i.valuation(K.factor_integer(2)[0][0])
            5
            sage: i.valuation(K.factor_integer(3)[0][0])
            0
            sage: i.valuation(0)
            Traceback (most recent call last):
            ...
            ValueError: p (= Fractional ideal (0) of Number Field in a with defining polynomial x^5 + 2) must be a nonzero prime
        """
        if not isinstance(p, NumberFieldIdeal):
            p = self.number_field().ideal(p)
        if p.is_zero() or not p.is_prime():
            raise ValueError, "p (= %s) must be a nonzero prime"%p
        if p.ring() != self.number_field():
            raise ValueError, "p (= %s) must be an ideal in %s"%self.number_field()
        nf = self.number_field().pari_nf()
        return ZZ(nf.idealval(self.pari_hnf(), p.__pari_prime))



def is_pari_zero_vector(z):
    """
    Return True if each entry of the PARI matrix row or vector z is 0.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field_ideal import is_pari_zero_vector
        sage: is_pari_zero_vector(pari('[]~'))
        True
        sage: is_pari_zero_vector(pari('[0,0]~'))
        True
        sage: is_pari_zero_vector(pari('[0,0,0,0,0]~'))
        True
        sage: is_pari_zero_vector(pari('[0,0,0,1,0]~'))
        False
    """
    for a in z:
        if a:
            return False
    return True
