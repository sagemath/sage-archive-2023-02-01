"""
Number Field Ideals

AUTHORS:

- Steven Sivek (2005-05-16)

- William Stein (2007-09-06): vastly improved the doctesting

- William Stein and John Cremona (2007-01-28): new class
  NumberFieldFractionalIdeal now used for all except the 0 ideal

We test that pickling works::

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

import sage.misc.latex as latex

import sage.rings.field_element as field_element
import sage.rings.polynomial.polynomial_element as polynomial
import sage.rings.polynomial.polynomial_ring as polynomial_ring
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.rational as rational
import sage.rings.integer as integer
import sage.rings.arith as arith
import sage.misc.misc as misc
from sage.rings.finite_field import FiniteField

import number_field
import number_field_element

from sage.libs.all import pari_gen
from sage.rings.ideal import (Ideal_generic, Ideal_fractional)
from sage.misc.misc import prod
from sage.structure.element import generic_power
from sage.structure.factorization import Factorization
from sage.structure.sequence import Sequence

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()

def convert_from_zk_basis(field, hnf):
    """
    Used internally in the number field ideal implementation for
    converting from the order basis to the number field basis.

    INPUT:

    -  ``field`` - a number field

    -   ``hnf`` - a pari HNF matrix, output by the pari_hnf() function.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_ideal import convert_from_zk_basis
        sage: k.<a> = NumberField(x^2 + 23)
        sage: I = k.factor(3)[0][0]; I
        Fractional ideal (3, -1/2*a + 1/2)
        sage: hnf = I.pari_hnf(); hnf
        [3, 0; 0, 1]
        sage: convert_from_zk_basis(k, hnf)
        [3, 1/2*x - 1/2]

    """
    return field.pari_nf().getattr('zk') * hnf

class NumberFieldIdeal(Ideal_generic):
    """
    An ideal of a number field.
    """
    def __init__(self, field, gens, coerce=True):
        """
        INPUT:

        -  ``field`` - a number field

        -   ``x`` - a list of NumberFieldElements belonging to the field

        EXAMPLES::

            sage: NumberField(x^2 + 1, 'a').ideal(7)
            Fractional ideal (7)
        """
        if not isinstance(field, number_field.NumberField_generic):
            raise TypeError, "field (=%s) must be a number field."%field

        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        if len(gens)==0:
            raise ValueError, "gens must have length at least 1 (zero ideal is not a fractional ideal)"
        Ideal_generic.__init__(self, field, gens, coerce)

    def _latex_(self):
        """
        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: K.ideal([2, 1/2*a - 1/2])._latex_()
            '\\left(2, \\frac{1}{2} a - \\frac{1}{2}\\right)'
            sage: latex(K.ideal([2, 1/2*a - 1/2]))
            \left(2, \frac{1}{2} a - \frac{1}{2}\right)
        """
        return '\\left(%s\\right)'%(", ".join([latex.latex(g) for g in \
                                                 self.gens_reduced()]))

    def __cmp__(self, other):
        """
        Compare an ideal of a number field to something else.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 3); K
            Number Field in a with defining polynomial x^2 + 3
            sage: f = K.factor(15); f
            (Fractional ideal (1/2*a - 3/2))^2 * (Fractional ideal (5))
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

    def coordinates(self, x):
        r"""
        Returns the coordinate vector of `x` with respect to this ideal.

        INPUT:
            ``x`` -- an element of the number field (or ring of integers) of this ideal.

        OUTPUT:
            A vector of length `n` (the degree of the field) giving
            the coordinates of `x` with respect to the integral basis
            of the ideal.  In general this will be a vector of
            rationals; it will consist of integers if and only if `x`
            is in the ideal.

        AUTHOR: John Cremona  2008-10-31
            Uses linear algebra.  The change-of-basis matrix is
            cached.  Provides simpler implementations for
            _contains_(), is_integral() and smallest_integer().

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: I = K.ideal(7+3*i)
            sage: Ibasis = I.integral_basis(); Ibasis
            [58, i + 41]
            sage: a = 23-14*i
            sage: acoords = I.coordinates(a); acoords
            (597/58, -14)
            sage: sum([Ibasis[j]*acoords[j] for j in range(2)]) == a
            True
            sage: b = 123+456*i
            sage: bcoords = I.coordinates(b); bcoords
            (-18573/58, 456)
            sage: sum([Ibasis[j]*bcoords[j] for j in range(2)]) == b
            True
       """
        K = self.number_field()
        V, from_V, to_V = K.absolute_vector_space()

        try:
            M = self.__basis_matrix_inverse
        except AttributeError:
            from sage.matrix.constructor import Matrix
            self.__basis_matrix_inverse = Matrix([to_V(b) for b in self.integral_basis()]).inverse()
            M = self.__basis_matrix_inverse
        return to_V(K(x))*M

    def _contains_(self, x):
        """
        Return True if x is an element of this ideal.

        This function is called (indirectly) when the \code{in}
        operator is used.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23); K
            Number Field in a with defining polynomial x^2 + 23
            sage: I = K.factor(13)[0][0]; I
            Fractional ideal (13, a - 4)
            sage: I._contains_(a)
            False
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
            sage: I = K.factor(13)[0][0]
            sage: I  # random sign in output
            Fractional ideal (-2*a^2 - 1)
            sage: 2/3 in I
            False
            sage: 1 in I
            False
            sage: 13 in I
            True
            sage: 1 in I*I^(-1)
            True
            sage: I   # random sign in output
            Fractional ideal (-2*a^2 - 1)
        """
        return self.coordinates(self.number_field()(x)).denominator() == 1

    def __elements_from_hnf(self, hnf):
        """
        Convert a PARI Hermite normal form matrix to a list of
        NumberFieldElements.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 389); K
            Number Field in a with defining polynomial x^3 + 389
            sage: I = K.factor(17)[0][0]
            sage: I       # random sign in generator
            Fractional ideal (-100*a^2 + 730*a - 5329)
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

    def __repr__(self):
        """
        Return the string representation of this number field ideal.

        NOTE: Only the zero ideal actually has type NumberFieldIdeal;
        all others have type NumberFieldFractionalIdeal.  So this
        function will only ever be called on the zero ideal.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: I = K.ideal(0); I
            Ideal (0) of Number Field in a with defining polynomial x^3 - 2
            sage: type(I)
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldIdeal'>
            sage: I = K.ideal(1); I
            Fractional ideal (1)
            sage: type(I)
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal'>
            sage: I = K.ideal(a); I
            Fractional ideal (a)
            sage: type(I)
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal'>
            sage: I = K.ideal(1/a); I
            Fractional ideal (1/2*a^2)
            sage: type(I)
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal'>
        """
        return "Ideal %s of %s"%(self._repr_short(), self.number_field())

    def _repr_short(self):
        """
        Efficient string representation of this fraction ideal.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 + 389); K
            Number Field in a with defining polynomial x^4 + 389
            sage: I = K.factor(17)[0][0]; I
            Fractional ideal (17, a^2 - 6)
            sage: I._repr_short()
            '(17, a^2 - 6)'
        """
        #NOTE -- we will *have* to not reduce the gens soon, since this
        # makes things insanely slow in general.
        # When I fix this, I *have* to also change the _latex_ method.
        return '(%s)'%(', '.join([str(x) for x in self.gens_reduced()]))
#        return '(%s)'%(', '.join([str(x) for x in self.gens()]))

    def _pari_(self):
        """
        Returns PARI Hermite Normal Form representations of this
        ideal.

        EXAMPLES::

            sage: K.<w> = NumberField(x^2 + 23)
            sage: I = K.class_group().0.ideal(); I
            Fractional ideal (2, 1/2*w - 1/2)
            sage: I._pari_()
            [2, 0; 0, 1]
        """
        return self.pari_hnf()

    def _pari_init_(self):
        """
        Returns self in PARI Hermite Normal Form as a string

        EXAMPLES::

            sage: K.<w> = NumberField(x^2 + 23)
            sage: I = K.class_group().0.ideal()
            sage: I._pari_init_()
            '[2, 0; 0, 1]'
        """
        return str(self._pari_())

    def pari_hnf(self):
        """
        Return PARI's representation of this ideal in Hermite normal form.

        EXAMPLES::

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

    def basis(self):
        """
        Return an immutable sequence of elements of this ideal (note:
        their parent is the number field) that form a basis for this
        ideal viewed as a ZZ-module.

        OUTPUT:
            basis -- an immutable sequence.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(7)
            sage: I = K.factor(11)[0][0]
            sage: I.basis()           # warning -- choice of basis can be somewhat random
            [11, 11*z, 11*z^2, z^3 + 5*z^2 + 4*z + 10, z^4 + z^2 + z + 5, z^5 + z^4 + z^3 + 2*z^2 + 6*z + 5]

        An example of a non-integral ideal.::

            sage: J = 1/I
            sage: J          # warning -- choice of generators can be somewhat random
            Fractional ideal (2/11*z^5 + 2/11*z^4 + 3/11*z^3 + 2/11)
            sage: J.basis()           # warning -- choice of basis can be somewhat random
            [1, z, z^2, 1/11*z^3 + 7/11*z^2 + 6/11*z + 10/11, 1/11*z^4 + 1/11*z^2 + 1/11*z + 7/11, 1/11*z^5 + 1/11*z^4 + 1/11*z^3 + 2/11*z^2 + 8/11*z + 7/11]
        """
        try:
            return self.__basis
        except AttributeError:
            pass
        hnf = self.pari_hnf()
        v = self.__elements_from_hnf(hnf)
        O = self.number_field().maximal_order()
        self.__basis = Sequence(v, immutable=True)
        return self.__basis

    def free_module(self):
        """
        Return the free ZZ-module contained in the vector space
        associated to the ambient number field, that corresponds
        to this ideal.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(7)
            sage: I = K.factor(11)[0][0]; I
            Fractional ideal (-2*z^4 - 2*z^2 - 2*z + 1)
            sage: A = I.free_module()
            sage: A              # warning -- choice of basis can be somewhat random
            Free module of degree 6 and rank 6 over Integer Ring
            User basis matrix:
            [11  0  0  0  0  0]
            [ 0 11  0  0  0  0]
            [ 0  0 11  0  0  0]
            [10  4  5  1  0  0]
            [ 5  1  1  0  1  0]
            [ 5  6  2  1  1  1]

        However, the actual ZZ-module is not at all random::

            sage: A.basis_matrix().change_ring(ZZ).echelon_form()
            [ 1  0  0  5  1  1]
            [ 0  1  0  1  1  7]
            [ 0  0  1  7  6 10]
            [ 0  0  0 11  0  0]
            [ 0  0  0  0 11  0]
            [ 0  0  0  0  0 11]

        The ideal doesn't have to be integral::

            sage: J = I^(-1)
            sage: B = J.free_module()
            sage: B.echelonized_basis_matrix()
            [ 1/11     0     0  7/11  1/11  1/11]
            [    0  1/11     0  1/11  1/11  5/11]
            [    0     0  1/11  5/11  4/11 10/11]
            [    0     0     0     1     0     0]
            [    0     0     0     0     1     0]
            [    0     0     0     0     0     1]

        This also works for relative extensions::

            sage: K.<a,b> = NumberField([x^2 + 1, x^2 + 2])
            sage: I = K.fractional_ideal(4)
            sage: I.free_module()
            Free module of degree 4 and rank 4 over Integer Ring
            User basis matrix:
            [  4   0   0   0]
            [ -3   7  -1   1]
            [  3   7   1   1]
            [  0 -10   0  -2]
            sage: J = I^(-1); J.free_module()
            Free module of degree 4 and rank 4 over Integer Ring
            User basis matrix:
            [  1/4     0     0     0]
            [-3/16  7/16 -1/16  1/16]
            [ 3/16  7/16  1/16  1/16]
            [    0  -5/8     0  -1/8]

        An example of intersecting ideals by intersecting free modules.::

            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: I = K.factor(2)
            sage: p1 = I[0][0]; p2 = I[1][0]
            sage: N = p1.free_module().intersection(p2.free_module()); N
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [  1 1/2 1/2]
            [  0   1   1]
            [  0   0   2]
            sage: N.index_in(p1.free_module()).abs()
            2

        TESTS:

        Sage can find the free module associated to quite large ideals
        quickly (see trac #4627)::

            sage: y = polygen(ZZ)
            sage: M.<a> = NumberField(y^20 - 2*y^19 + 10*y^17 - 15*y^16 + 40*y^14 - 64*y^13 + 46*y^12 + 8*y^11 - 32*y^10 + 8*y^9 + 46*y^8 - 64*y^7 + 40*y^6 - 15*y^4 + 10*y^3 - 2*y + 1)
            sage: M.ideal(prod(prime_range(6000, 6200))).free_module()
            Free module of degree 20 and rank 20 over Integer Ring
            User basis matrix:
            20 x 20 dense matrix over Rational Field

        """
        try:
            return self.__free_module
        except AttributeError:
            pass
        M = basis_to_module(self.basis(), self.number_field())
        self.__free_module = M
        return M

    def reduce_equiv(self):
        """
        Return a small ideal that is equivalent to self in the group
        of fractional ideals modulo principal ideals.  Very often (but
        not always) if self is principal then this function returns
        the unit ideal.

        ALGORITHM: Calls pari's idealred function.

        EXAMPLES::

            sage: K.<w> = NumberField(x^2 + 23)
            sage: I = ideal(w*23^5); I
            Fractional ideal (6436343*w)
            sage: I.reduce_equiv()
            Fractional ideal (1)
            sage: I = K.class_group().0.ideal()^10; I
            Fractional ideal (1024, 1/2*w + 979/2)
            sage: I.reduce_equiv()
            Fractional ideal (2, 1/2*w - 1/2)
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
        ``idealtwoelt`` function, which takes exponential time, the
        first time it is called for each ideal.  Also, this indirectly
        uses ``bnfisprincipal``, so set ``proof=True`` if you
        want to prove correctness (which *is* the default).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1, 'i')
            sage: J = K.ideal([i+1, 2])
            sage: J.gens()
            (i + 1, 2)
            sage: J.gens_reduced()
            (i + 1,)

        TESTS::

            sage: all(j.parent() is K for j in J.gens())
            True
            sage: all(j.parent() is K for j in J.gens_reduced())
            True

            sage: K.<a> = NumberField(x^4 + 10*x^2 + 20)
            sage: J = K.prime_above(5)
            sage: J.is_principal()
            False
            sage: J.gens_reduced()
            (5, a)
            sage: all(j.parent() is K for j in J.gens())
            True
            sage: all(j.parent() is K for j in J.gens_reduced())
            True
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
            gens = [ K(a), K(R(nf.getattr('zk')*alpha)) ]
            if gens[1] in K.ideal(gens[0]):
                gens = gens[:1]
            elif gens[0] in K.ideal(gens[1]):
                gens = gens[1:]
            self.__reduced_generators = tuple(gens)
            return self.__reduced_generators

    def integral_basis(self):
        r"""
        Return a list of generators for this ideal as a `\mathbb{Z}`-module.

        EXAMPLES::

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

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2-5)
            sage: I = K.ideal(2/(5+a))
            sage: I.is_integral()
            False
            sage: J,d = I.integral_split()
            sage: J
            Fractional ideal (-1/2*a + 5/2)
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

        EXAMPLES::

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
            self.__is_integral = all([a in one for a in self.integral_basis()])
            return self.__is_integral

    def is_maximal(self):
        """
        Return True if this ideal is maximal.  This is equivalent to
        self being prime and nonzero.

        EXAMPLES::

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

        EXAMPLES::

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
            return self._pari_prime is not None
        except AttributeError:
            K = self.number_field()
            F = list(K.pari_nf().idealfactor(self.pari_hnf()))
            ### We should definitely cache F as the factorization of self
            P, exps = F[0], F[1]
            if len(P) != 1 or exps[0] != 1:
                self._pari_prime = None
            else:
                self._pari_prime = P[0]
            return self._pari_prime is not None

    def is_principal(self, proof=None):
        r"""
        Return True if this ideal is principal.

        Since it uses the PARI method \code{bnfisprincipal}, specify
        \code{proof=True} (this is the default setting) to prove the
        correctness of the output.

        EXAMPLES:

        We create equal ideals in two different ways, and note that
        they are both actually principal ideals.::

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
            self.__is_principal = not any(v[0])
            if self.__is_principal:
                K = self.number_field()
                R = K.polynomial().parent()
                g = K(R(bnf.getattr('zk') * v[1]))
                self.__reduced_generators = tuple([g])
            return self.__is_principal

    def is_zero(self):
        """
        Return True iff self is the zero ideal

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 2); K
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(3).is_zero()
            False
            sage: I=K.ideal(0); I.is_zero()
            True
            sage: I
            Ideal (0) of Number Field in a with defining polynomial x^2 + 2

            (0 is a NumberFieldIdeal, not a NumberFieldFractionIdeal)
        """
        return self == self.number_field().ideal(0)

    def norm(self):
        """
        Return the norm of this fractional ideal as a rational number.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 + 23); K
            Number Field in a with defining polynomial x^4 + 23
            sage: I = K.ideal(19); I
            Fractional ideal (19)
            sage: factor(I.norm())
            19^4
            sage: F = I.factor()
            sage: F[0][0].norm().factor()
            19^2
        """
        try:
            return self._norm
        except AttributeError:
            pass
        self._norm = QQ(self.number_field().pari_nf().idealnorm(self.pari_hnf()))
        return self._norm

    # synonyms (using terminology of relative number fields)

    def absolute_norm(self):
        """
        A synonym for norm.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.ideal(1 + 2*i).absolute_norm()
            5
        """
        return self.norm()

    def relative_norm(self):
        """
        A synonym for norm.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.ideal(1 + 2*i).relative_norm()
            5
        """
        return self.norm()

    def absolute_ramification_index(self):
        """
        A synonym for ramification_index.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.ideal(1 + i).absolute_ramification_index()
            2
        """
        return self.ramification_index()

    def relative_ramification_index(self):
        """
        A synonym for ramification_index.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.ideal(1 + i).relative_ramification_index()
            2
        """
        return self.ramification_index()

    def number_field(self):
        """
        Return the number field that this is a fractional ideal in.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 2); K
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(3).number_field()
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(0).number_field() # not tested (not implemented)
            Number Field in a with defining polynomial x^2 + 2
        """
        return self.ring()

    def smallest_integer(self):
        r"""
        Return the smallest non-negative integer in `I \cap \mathbb{Z}`,
        where `I` is this ideal.  If `I = 0`, returns 0.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2+6)
            sage: I = K.ideal([4,a])/7
            sage: I.smallest_integer()
            2

            sage: K.<i> = QuadraticField(-1)
            sage: P1, P2 = [P for P,e in K.factor(13)]
            sage: all([(P1^i*P2^j).smallest_integer() == 13^max(i,j,0) for i in range(-3,3) for j in range(-3,3)])
            True
            sage: I = K.ideal(0)
            sage: I.smallest_integer()
            0

            # See trac\# 4392:
            sage: K.<a>=QuadraticField(-5)
            sage: I=K.ideal(7)
            sage: I.smallest_integer()
            7

            sage: K.<z> = CyclotomicField(13)
            sage: a = K([-8, -4, -4, -6, 3, -4, 8, 0, 7, 4, 1, 2])
            sage: I = K.ideal(a)
            sage: I.smallest_integer()
            146196692151
            sage: I.norm()
            1315770229359
            sage: I.norm() / I.smallest_integer()
            9
        """
        try:
            return self.__smallest_integer
        except AttributeError:
            pass

        if self.is_zero():
            self.__smallest_integer = ZZ(0)
            return self.__smallest_integer

        if self.is_prime():
            self.__smallest_integer = ZZ(self._pari_prime.getattr('p'))
            return self.__smallest_integer


#   New code by John Cremona, 2008-10-30, using the new coordinates()
#   function instead of factorization.
#
#   Idea: We write 1 as a Q-linear combination of the Z-basis of self,
#         and return the denominator of this vector.
#
#   Note: It seems that in practice for integral ideals the first
#         element of the integral basis is the smallest integer, but
#         we cannot rely on this.

        self.__smallest_integer =  self.coordinates(1).denominator()
        return self.__smallest_integer

    def valuation(self, p):
        r"""
        Return the valuation of this fractional ideal at the prime
        `\mathfrak{p}`.  If `\mathfrak{p}` is not prime, raise a
        ValueError.

        INPUT:
            p -- a prime ideal of this number field.

        OUTPUT:
            integer

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 + 2); K
            Number Field in a with defining polynomial x^5 + 2
            sage: i = K.ideal(38); i
            Fractional ideal (38)
            sage: i.valuation(K.factor(19)[0][0])
            1
            sage: i.valuation(K.factor(2)[0][0])
            5
            sage: i.valuation(K.factor(3)[0][0])
            0
            sage: i.valuation(0)
            Traceback (most recent call last):
            ...
            ValueError: p (= 0) must be nonzero
        """
        if p==0:
            raise ValueError, "p (= %s) must be nonzero"%p
        if not isinstance(p, NumberFieldFractionalIdeal):
            p = self.number_field().ideal(p)
        if not p.is_prime():
            raise ValueError, "p (= %s) must be a prime"%p
        if p.ring() != self.number_field():
            raise ValueError, "p (= %s) must be an ideal in %s"%self.number_field()
        nf = self.number_field().pari_nf()
        return ZZ(nf.idealval(self.pari_hnf(), p._pari_prime))

    def decomposition_group(self):
        """
        Return the decomposition group of self, as a subset of the automorphism
        group of the number field of self. Raises an error if the field isn't
        Galois. See the decomposition_group method of the GaloisGroup_v2 class
        for further examples and doctests.

        EXAMPLE::

            sage: QuadraticField(-23, 'w').primes_above(7)[0].decomposition_group()
            Galois group of Number Field in w with defining polynomial x^2 + 23
        """
        return self.number_field().galois_group().decomposition_group(self)

    def ramification_group(self, v):
        """
        Return the vth ramification group of self, i.e. the set of elements s
        of the Galois group of the number field of self (which we assume is Galois)
        such that s acts trivially modulo self^(v+1). See the
        ramification_group method of the GaloisGroup class for further examples
        and doctests.

        EXAMPLE::

            sage: QuadraticField(-23, 'w').primes_above(23)[0].ramification_group(0)
            Galois group of Number Field in w with defining polynomial x^2 + 23
            sage: QuadraticField(-23, 'w').primes_above(23)[0].ramification_group(1)
            Subgroup [()] of Galois group of Number Field in w with defining polynomial x^2 + 23
        """

        return self.number_field().galois_group().ramification_group(self, v)

    def inertia_group(self):
        """
        Return the inertia group of self, i.e. the set of elements s of the
        Galois group of the number field of self (which we assume is Galois)
        such that s acts trivially modulo self. This is the same as the 0th
        ramification group of self. See the inertia_group method of the
        GaloisGroup_v2 class for further examples and doctests.

        EXAMPLE::

            sage: QuadraticField(-23, 'w').primes_above(23)[0].inertia_group()
            Galois group of Number Field in w with defining polynomial x^2 + 23
        """
        return self.ramification_group(0)

    def artin_symbol(self):
        """
        Return the Artin symbol ( K / Q, self), where K is the number field of self.
        This is the unique element s of the decomposition group of self such that
        s(x) = x^p mod self where p is the residue characteristic of self.
        (Here self should be prime and unramified.)

        See the artin_symbol method of the GaloisGroup_v2 class for further
        documentation and examples.

        EXAMPLE::

            sage: QuadraticField(-23, 'w').primes_above(7)[0].artin_symbol()
            (1,2)
        """
        return self.number_field().galois_group().artin_symbol(self)

def basis_to_module(B, K):
    """
    Given a basis B of elements for a ZZ-submodule of a number field K, return
    the corresponding ZZ-submodule.

    EXAMPLES::

        sage: K.<w> = NumberField(x^4 + 1)
        sage: from sage.rings.number_field.number_field_ideal import basis_to_module
        sage: basis_to_module([K.0, K.0^2 + 3], K)
        Free module of degree 4 and rank 2 over Integer Ring
        User basis matrix:
        [0 1 0 0]
        [3 0 1 0]
    """
    V, from_V, to_V = K.absolute_vector_space()
    M = ZZ**(V.dimension())
    C = [to_V(K(b)) for b in B]
    return M.span_of_basis(C)

def is_NumberFieldIdeal(x):
    """
    Return True if x is an ideal of a number field.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_ideal import is_NumberFieldIdeal
        sage: is_NumberFieldIdeal(2/3)
        False
        sage: is_NumberFieldIdeal(ideal(5))
        False
        sage: k.<a> = NumberField(x^2 + 2)
        sage: I = k.ideal([a + 1]); I
        Fractional ideal (a + 1)
        sage: is_NumberFieldIdeal(I)
        True
        sage: Z = k.ideal(0); Z
        Ideal (0) of Number Field in a with defining polynomial x^2 + 2
        sage: is_NumberFieldIdeal(Z)
        True
    """
    return isinstance(x, NumberFieldIdeal)


class NumberFieldFractionalIdeal(NumberFieldIdeal):

    def __init__(self, field, gens, coerce=True):
        """
        INPUT:
            field -- a number field
            x -- a list of NumberFieldElements of the field, not all zero

        EXAMPLES::

            sage: NumberField(x^2 + 1, 'a').ideal(7)
            Fractional ideal (7)
        """
        if not isinstance(field, number_field.NumberField_generic):
            raise TypeError, "field (=%s) must be a number field."%field

        if len(gens)==0:
            raise ValueError, "gens must have length at least 1 (zero ideal is not a fractional ideal)"
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        if misc.exists(gens,bool)[0]:
            NumberFieldIdeal.__init__(self, field, gens)
        else:
            raise ValueError, "gens must have a nonzero element (zero ideal is not a fractional ideal)"

    def __repr__(self):
        """
        Return the string representation of this number field fractional ideal.

        NOTE: Only the zero ideal actually has type NumberFieldIdeal;
        all others have type NumberFieldFractionalIdeal.

        EXAMPLES::

            sage: K.<a>=NumberField(x^2+5)
            sage: I = K.ideal([2,1+a]); I
            Fractional ideal (2, a + 1)
            sage: type(I)
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal'>
        """
        return "Fractional ideal %s"%self._repr_short()

    def divides(self, other):
        """
        Returns True if this ideal divides other and False otherwise.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(11); K
            Cyclotomic Field of order 11 and degree 10
            sage: I = K.factor(31)[0][0]; I
            Fractional ideal (-3*a^7 - 4*a^5 - 3*a^4 - 3*a^2 - 3*a - 3)
            sage: I.divides(I)
            True
            sage: I.divides(31)
            True
            sage: I.divides(29)
            False
        """
        if not isinstance(other, NumberFieldIdeal):
            other = self.number_field().ideal(other)
        return (other / self).is_integral()

    def factor(self):
        """
        Factorization of this ideal in terms of prime ideals.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 + 23); K
            Number Field in a with defining polynomial x^4 + 23
            sage: I = K.ideal(19); I
            Fractional ideal (19)
            sage: F = I.factor(); F
            (Fractional ideal (a^2 + 2*a + 2)) * (Fractional ideal (a^2 - 2*a + 2))
            sage: type(F)
            <class 'sage.structure.factorization.Factorization'>
            sage: list(F)
            [(Fractional ideal (a^2 + 2*a + 2), 1), (Fractional ideal (a^2 - 2*a + 2), 1)]
            sage: F.prod()
            Fractional ideal (19)
        """
        try:
            return self.__factorization
        except AttributeError:
            K = self.number_field()
            F = list(K.pari_nf().idealfactor(self.pari_hnf()))
            P, exps = F[0], F[1]
            A = []
            zk_basis = K.pari_nf().getattr('zk')
            for i, p in enumerate(P):
                prime, alpha = p.getattr('gen')
                I = K.ideal([ZZ(prime), K(zk_basis * alpha)])
                I._pari_prime = p
                A.append((I,ZZ(exps[i])))
            self.__factorization = Factorization(A)
            return self.__factorization

    def prime_factors(self):
        """
        Return a list of the prime ideal factors of self

        OUTPUT:
            list -- list of prime ideals (a new list is returned
            each time this function is called)

        EXAMPLES::

            sage: K.<w> = NumberField(x^2 + 23)
            sage: I = ideal(w+1)
            sage: I.prime_factors()
            [Fractional ideal (2, 1/2*w - 1/2),
            Fractional ideal (2, 1/2*w + 1/2),
            Fractional ideal (3, -1/2*w - 1/2)]
        """
        return [x[0] for x in self.factor()]

    def __div__(self, other):
        """
        Return the quotient self / other.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2 - 5)
            sage: I = K.ideal(2/(5+a))
            sage: J = K.ideal(17+a)
            sage: I/J
            Fractional ideal (-17/1420*a + 1/284)
            sage: (I/J) * J
            Fractional ideal (-1/5*a)
            sage: (I/J) * J == I
            True
        """
        return self * other.__invert__()

    def __invert__(self):
        """
        Return the multiplicative inverse of self.  Call with ~self.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: I = K.ideal(2/(5+a))
            sage: ~I
            Fractional ideal (1/2*a + 5/2)
            sage: 1/I
            Fractional ideal (1/2*a + 5/2)
            sage: (1/I) * I
            Fractional ideal (1)
        """
        nf = self.number_field().pari_nf()
        hnf = nf.idealdiv(self.number_field().ideal(1).pari_hnf(),
                          self.pari_hnf())
        I = self.number_field().ideal(NumberFieldIdeal._NumberFieldIdeal__elements_from_hnf(self,hnf))
        I.__pari_hnf = hnf
        return I

    def __pow__(self, r):
        """
        Return self to the power of r.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: I = K.ideal(2/(5+a))
            sage: J = I^2
            sage: Jinv = I^(-2)
            sage: J*Jinv
            Fractional ideal (1)
        """
        return generic_power(self, r)

    def is_maximal(self):
        """
        Return True if this ideal is maximal.  This is equivalent to
        self being prime, since it is nonzero.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 3); K
            Number Field in a with defining polynomial x^3 + 3
            sage: K.ideal(5).is_maximal()
            False
            sage: K.ideal(7).is_maximal()
            True
        """
        return self.is_prime()

    def is_trivial(self, proof=None):
        """
        Returns True if this is a trivial ideal.

        EXAMPLES::

            sage: F.<a> = QuadraticField(-5)
            sage: I = F.ideal(3)
            sage: I.is_trivial()
            False
            sage: J = F.ideal(5)
            sage: J.is_trivial()
            False
            sage: (I+J).is_trivial()
            True
        """
        return self == self.number_field().ideal(1)

    def ramification_index(self):
        r"""
        Return the ramification index of this fractional ideal,
        assuming it is prime.  Otherwise, raise a ValueError.

        The ramification index is the power of this prime appearing in
        the factorization of the prime in `\mathbb{Z}` that this prime lies
        over.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 2); K
            Number Field in a with defining polynomial x^2 + 2
            sage: f = K.factor(2); f
            (Fractional ideal (-a))^2
            sage: f[0][0].ramification_index()
            2
            sage: K.ideal(13).ramification_index()
            1
            sage: K.ideal(17).ramification_index()
            Traceback (most recent call last):
            ...
            ValueError: the fractional ideal (= Fractional ideal (17)) is not prime
        """
        if self.is_prime():
            return ZZ(self._pari_prime.getattr('e'))
        raise ValueError, "the fractional ideal (= %s) is not prime"%self

    def residues(self):
        """
        Returns a iterator through a complete list of residues modulo this integral ideal.

        An error is raised if this fractional ideal is not integral.

        AUTHOR: John Cremona

        EXAMPLES::


            sage: K.<i>=NumberField(x^2+1)
            sage: res =  K.ideal(2).residues(); res  # random address
            xmrange_iter([[0, 1], [0, 1]], <function <lambda> at 0xa252144>)
            sage: list(res)
            [0, 1, i, i + 1]
            sage: list(K.ideal(2+i).residues())
            [-2, -1, 0, 1, 2]
            sage: list(K.ideal(i).residues())
            [0]
            sage: I = K.ideal(3+6*i)
            sage: reps=I.residues()
            sage: len(list(reps)) == I.norm()
            True
            sage: all([r==s or not (r-s) in I for r in reps for s in reps])
            True

            sage: K.<a> = NumberField(x^3-10)
            sage: I = K.ideal(a-1)
            sage: len(list(I.residues())) == I.norm()
            True

            sage: K.<z> = CyclotomicField(11)
            sage: len(list(K.primes_above(3)[0].residues())) == 3**5
            True
        """
        R = self.number_field().maximal_order()
        Rbasis = R.basis()
        n = len(Rbasis)
        from sage.matrix.all import MatrixSpace
        try:
            M = MatrixSpace(ZZ,n)([R.coordinates(y) for y in self.basis()])
        except TypeError:
            raise TypeError, "residues only defined for integral ideals"

        D, U, V = M.smith_form()
        from sage.misc.mrange import xmrange_iter
        d = [D[i,i] for i in range(n)]
        coord_ranges = [range((-di+2)//2,(di+2)//2) for di in d]
        V = V.inverse()
        newRbasis = [sum([V[i,j]*Rbasis[j] for j in range(n)]) for i in range(n)]
        combo = lambda c: sum([c[i]*newRbasis[i] for i in range(n)])
        return xmrange_iter(coord_ranges, combo)

    def invertible_residues(self):
        """
        Returns a iterator through a list of invertible residues modulo this integral ideal.

        An error is raised if this fractional ideal is not integral.

        AUTHOR: John Cremona

        EXAMPLES::

            sage: K.<i>=NumberField(x^2+1)
            sage: ires =  K.ideal(2).invertible_residues(); ires  # random address
            <generator object at 0xa2feb6c>
            sage: list(ires)
            [1, i]
            sage: list(K.ideal(2+i).invertible_residues())
            [-2, -1, 1, 2]
            sage: list(K.ideal(i).residues())
            [0]
            sage: I = K.ideal(3+6*i)
            sage: units=I.invertible_residues()
            sage: len(list(units))==I.euler_phi()
            True

            sage: K.<a> = NumberField(x^3-10)
            sage: I = K.ideal(a-1)
            sage: len(list(I.invertible_residues())) == I.euler_phi()
            True

            sage: K.<z> = CyclotomicField(10)
            sage: len(list(K.primes_above(3)[0].invertible_residues()))
            80
        """
        for r in self.residues():
            if self.is_coprime(r):
                yield r


    def denominator(self):
        r"""
        Return the denominator ideal of this fractional ideal. Each fractional
        ideal has a unique expression as `N/D` where N, D are coprime integral
        ideals; the denominator is D.

        EXAMPLES::

            sage: K.<i>=NumberField(x^2+1)
            sage: I = K.ideal((3+4*i)/5); I
            Fractional ideal (4/5*i + 3/5)
            sage: I.denominator()
            Fractional ideal (2*i + 1)
            sage: I.numerator()
            Fractional ideal (-i - 2)
            sage: I.numerator().is_integral() and I.denominator().is_integral()
            True
            sage: I.numerator() + I.denominator() == K.unit_ideal()
            True
            sage: I.numerator()/I.denominator() == I
            True
        """
        try:
            return self._denom_ideal
        except AttributeError:
            pass
        self._denom_ideal = (self + self.number_field().unit_ideal())**(-1)
        return self._denom_ideal

    def numerator(self):
        """
        Return the numerator ideal of this fractional ideal.

        Each fractional ideal has a unique expression as N/D where N,
        D are coprime integral ideals.  The numerator is N.

        EXAMPLES::

            sage: K.<i>=NumberField(x^2+1)
            sage: I = K.ideal((3+4*i)/5); I
            Fractional ideal (4/5*i + 3/5)
            sage: I.denominator()
            Fractional ideal (2*i + 1)
            sage: I.numerator()
            Fractional ideal (-i - 2)
            sage: I.numerator().is_integral() and I.denominator().is_integral()
            True
            sage: I.numerator() + I.denominator() == K.unit_ideal()
            True
            sage: I.numerator()/I.denominator() == I
            True
        """
        try:
            return self._num_ideal
        except AttributeError:
            pass
        self._num_ideal = self * self.denominator()
        return self._num_ideal

    def is_coprime(self, other):
        """
        Returns True if this ideal is coprime to the other, else False.

        INPUT:
           other -- another ideal of the same field, or generators of an ideal.

        OUTPUT:
           True if self and other are coprime, else False.

        NOTE:
           This function works for fractional ideals as well as
           integral ideals.

        AUTHOR: John Cremona

        EXAMPLES::

            sage: K.<i>=NumberField(x^2+1)
            sage: I = K.ideal(2+i)
            sage: J = K.ideal(2-i)
            sage: I.is_coprime(J)
            True
            sage: (I^-1).is_coprime(J^3)
            True
            sage: I.is_coprime(5)
            False
            sage: I.is_coprime(6+i)
            True

            # See trac \# 4536:
            sage: E.<a> = NumberField(x^5 + 7*x^4 + 18*x^2 + x - 3)
            sage: OE = E.ring_of_integers()
            sage: i,j,k = [u[0] for u in factor(3*OE)]
            sage: (i/j).is_coprime(j/k)
            False
            sage: (j/k).is_coprime(j/k)
            False

            sage: F.<a, b> = NumberField([x^2 - 2, x^2 - 3])
            sage: F.ideal(3 - a*b).is_coprime(F.ideal(3))
            False
        """
        # Catch invalid inputs by making sure that we can make an ideal out of other.
        K = self.number_field()
        one = K.unit_ideal()
        other = K.ideal(other)
        if self.is_integral() and other.is_integral():
            if arith.gcd(ZZ(self.absolute_norm()), ZZ(other.absolute_norm())) == 1:
                return True
            else:
                return self+other == one
        # This special case is necessary since the zero ideal is not a
        # fractional ideal!
        if other.absolute_norm() == 0:
            return self == one
        D1 = self.denominator()
        N1 = self.numerator()
        D2 = other.denominator()
        N2 = other.numerator()
        return N1+N2==one and N1+D2==one and D1+N2==one and D1+D2==one

    def element_1_mod(self, other):
        """
        Returns an element r in this ideal such that 1-r is in other

        An error is raised if either ideal is not integral of if they
        are not coprime.

        INPUT:
            other -- another ideal of the same field, or generators of an ideal.
        OUTPUT:
            r -- an element of the ideal self such that 1-r is in the ideal other

        AUTHOR: Maite Aranes

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: A = K.ideal(a+1); A; A.norm()
            Fractional ideal (a + 1)
            3
            sage: B = K.ideal(a^2-4*a+2); B; B.norm()
            Fractional ideal (a^2 - 4*a + 2)
            68
            sage: r = A.element_1_mod(B); r
            a^2 + 5
            sage: r in A
            True
            sage: 1-r in B
            True

        TESTS::

            sage: K.<a> = NumberField(x^3-2)
            sage: A = K.ideal(a+1)
            sage: B = K.ideal(a^2-4*a+1); B; B.norm()
            Fractional ideal (a^2 - 4*a + 1)
            99
            sage: A.element_1_mod(B)
            Traceback (most recent call last):
            ...
            TypeError: Fractional ideal (a + 1), Fractional ideal (a^2 - 4*a + 1) are not coprime ideals

            sage: B = K.ideal(1/a); B
            Fractional ideal (1/2*a^2)
            sage: A.element_1_mod(B)
            Traceback (most recent call last):
            TypeError: element_1_mod only defined for integral ideals
        """
        # Catch invalid inputs by making sure that we can make an ideal out of other.
        k = self.number_field()
        other = k.ideal(other)

        #we want a basis for the ring of integers with first element=1
        R = k.unit_ideal()
        Rbasis = R.basis()
        assert Rbasis[0]==1  # true in 3.2.2

        n = len(Rbasis)

        #matrices for self and other in terms of basis chosen for R
        self_b = self.basis()
        other_b = other.basis()

        from sage.matrix.all import MatrixSpace

        try:
            M_self = MatrixSpace(ZZ,n)([R.coordinates(y) for y in self_b])
            M_other = MatrixSpace(ZZ,n)([R.coordinates(y) for y in other_b])
        except TypeError:
            raise TypeError, "element_1_mod only defined for integral ideals"

        #hnf for matrix representing C = self+other:
        C = M_self.stack(M_other)
        Chnf, U = C.hermite_form(transformation=True)

        #we make sure the ideals self and other are coprime
        if Chnf[0][0]!=1 or not (Chnf.submatrix(0,1,1)).is_zero():
            raise TypeError, "%s, %s are not coprime ideals"%(self, other)

        #element r in self such that 1 - r in other
        from sage.modules.free_module_element import vector
        r = vector([U[0][i] for i in range(n)])*M_self
        r = sum([r[i]*Rbasis[i] for i in range(n)])

        return r


    def euler_phi(self):
        r"""
        Returns the Euler `\varphi`-function of this integral ideal.

        This is the order of the multiplicative group of the quotient
        modulo the ideal.

        An error is raised if the ideal is not integral.

        EXAMPLES::

            sage: K.<i>=NumberField(x^2+1)
            sage: I = K.ideal(2+i)
            sage: [r for r in I.residues() if I.is_coprime(r)]
            [-2, -1, 1, 2]
            sage: I.euler_phi()
            4
            sage: J = I^3
            sage: J.euler_phi()
            100
            sage: len([r for r in J.residues() if J.is_coprime(r)])
            100
            sage: J = K.ideal(3-2*i)
            sage: I.is_coprime(J)
            True
            sage: I.euler_phi()*J.euler_phi() == (I*J).euler_phi()
            True
            sage: L.<b> = K.extension(x^2 - 7)
            sage: L.ideal(3).euler_phi()
            64
        """
        if not self.is_integral():
            raise ValueError, "euler_phi only defined for integral ideals"
        return prod([(np-1)*np**(e-1) \
                     for np,e in [(p.absolute_norm(),e) \
                                  for p,e in self.factor()]])

    def prime_to_idealM_part(self, M):
        """
        Version for integral ideals of the prime_to_m_part function over ZZ.
        Returns the largest divisor of self that is coprime to the ideal M.

        INPUT:
            M -- an integral ideal of the same field, or generators of an ideal

        OUTPUT:
            An ideal which is the largest divisor of self that is coprime to M.

        AUTHOR: Maite Aranes

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 23)
            sage: I = k.ideal(a+1)
            sage: M = k.ideal(2, 1/2*a - 1/2)
            sage: J = I.prime_to_idealM_part(M); J
            Fractional ideal (12, 1/2*a + 13/2)
            sage: J.is_coprime(M)
            True

            sage: J = I.prime_to_idealM_part(2); J
            Fractional ideal (3, -1/2*a - 1/2)
            sage: J.is_coprime(M)
            True
        """
        # Catch invalid inputs by making sure that we can make an ideal out of M.
        k = self.number_field()
        M = k.ideal(M)

        if not self.is_integral or not M.is_integral():
            raise TypeError, "prime_to_idealM_part defined only for integral ideals"

        if self.is_coprime(M):
            return self
        G = self + M
        I = self
        while not G.is_trivial():
            I = I/G
            G = I + G
        return I

    def _p_quotient(self, p):
        """
        This is an internal technical function that is used for example for
        computing the quotient of the ring of integers by a prime ideal.

        INPUT:
            p -- a prime number contained in self.

        OUTPUT:
            V -- a vector space of characteristic p
            quo -- a partially defined quotient homomorphism from the
                   ambient number field to V
            lift -- a section of quo.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1); O = K.maximal_order()
            sage: I = K.factor(3)[0][0]
            sage: Q, quo, lift = I._p_quotient(3); Q
            Vector space quotient V/W of dimension 2 over Finite Field of size 3 where
            V: Vector space of dimension 2 over Finite Field of size 3
            W: Vector space of degree 2 and dimension 0 over Finite Field of size 3
            Basis matrix:
            []

        We do an example with a split prime and show both the quo and lift maps:
            sage: K.<i> = NumberField(x^2 + 1); O = K.maximal_order()
            sage: I = K.factor(5)[0][0]
            sage: Q,quo,lift = I._p_quotient(5)
            sage: lift(quo(i))
            3
            sage: lift(quo(i)) - i in I
            True
            sage: quo(lift(Q.0))
            (1)
            sage: Q.0
            (1)
            sage: Q
            Vector space quotient V/W of dimension 1 over Finite Field of size 5 where
            V: Vector space of dimension 2 over Finite Field of size 5
            W: Vector space of degree 2 and dimension 1 over Finite Field of size 5
            Basis matrix:
            [1 3]
            sage: quo
            Partially defined quotient map from Number Field in i with defining polynomial x^2 + 1 to an explicit vector space representation for the quotient of the ring of integers by (p,I) for the ideal I=Fractional ideal (-i - 2).
            sage: lift
            Lifting map to Maximal Order in Number Field in i with defining polynomial x^2 + 1 from quotient of integers by Fractional ideal (-i - 2)
        """
        return quotient_char_p(self, p)

    def residue_field(self, names=None):
        r"""
        Return the residue class field of this fractional ideal, which
        must be prime.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: P.residue_field()
            Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
            sage: P.residue_field('z')
            Residue field in z of Fractional ideal (2*a^2 + 3*a - 10)

        Another example::

            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(389).factor()[0][0]; P
            Fractional ideal (389, a^2 - 44*a - 9)
            sage: P.residue_class_degree()
            2
            sage: P.residue_field()
            Residue field in abar of Fractional ideal (389, a^2 - 44*a - 9)
            sage: P.residue_field('z')
            Residue field in z of Fractional ideal (389, a^2 - 44*a - 9)
            sage: FF.<w> = P.residue_field()
            sage: FF
            Residue field in w of Fractional ideal (389, a^2 - 44*a - 9)
            sage: FF((a+1)^390)
            36
            sage: FF(a)
            w

        An example of reduction maps to the residue field: these are
        defined on the whole valuation ring, i.e. the subring of the
        number field consisting of elements with non-negative
        valuation.  This shows that the issue raised in trac \#1951
        has been fixed.
            sage: K.<i> = NumberField(x^2 + 1)
            sage: P1, P2 = [g[0] for g in K.factor(5)]; (P1,P2)
            (Fractional ideal (-i - 2), Fractional ideal (2*i + 1))
            sage: a = 1/(1+2*i)
            sage: F1, F2 = [g.residue_field() for g in [P1,P2]]; (F1,F2)
            (Residue field of Fractional ideal (-i - 2),
            Residue field of Fractional ideal (2*i + 1))
            sage: a.valuation(P1)
            0
            sage: F1(i/7)
            4
            sage: F1(a)
            3
            sage: a.valuation(P2)
            -1
            sage: F2(a)
            Traceback (most recent call last):
            ZeroDivisionError: Cannot reduce field element -2/5*i + 1/5 modulo Fractional ideal (2*i + 1): it has negative valuation
        """
        if not self.is_prime():
            raise ValueError, "The ideal must be prime"
        return self.number_field().residue_field(self, names = names)

    def residue_class_degree(self):
        r"""
        Return the residue class degree of this fractional ideal,
        assuming it is prime.  Otherwise, raise a ValueError.

        The residue class degree of a prime ideal `I` is the degree of
        the extension `O_K/I` of its prime subfield.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 + 2); K
            Number Field in a with defining polynomial x^5 + 2
            sage: f = K.factor(19); f
            (Fractional ideal (a^2 + a - 3)) * (Fractional ideal (-2*a^4 - a^2 + 2*a - 1)) * (Fractional ideal (a^2 + a - 1))
            sage: [i.residue_class_degree() for i, _ in f]
            [2, 2, 1]
        """
        if self.is_prime():
            return ZZ(self._pari_prime.getattr('f'))
        raise ValueError, "the ideal (= %s) is not prime"%self

def is_NumberFieldFractionalIdeal(x):
    """
    Return True if x is a fractional ideal of a number field.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_ideal import is_NumberFieldFractionalIdeal
        sage: is_NumberFieldFractionalIdeal(2/3)
        False
        sage: is_NumberFieldFractionalIdeal(ideal(5))
        False
        sage: k.<a> = NumberField(x^2 + 2)
        sage: I = k.ideal([a + 1]); I
        Fractional ideal (a + 1)
        sage: is_NumberFieldFractionalIdeal(I)
        True
        sage: Z = k.ideal(0); Z
        Ideal (0) of Number Field in a with defining polynomial x^2 + 2
        sage: is_NumberFieldFractionalIdeal(Z)
        False
    """
    return isinstance(x, NumberFieldFractionalIdeal)

class QuotientMap:
    """
    Class to hold data needed by quotient maps from number field
    orders to residue fields.  These are only partial maps: the exact
    domain is the appropriate valuation ring.  For examples, see
    residue_field().
    """
    def __init__(self, K, M_OK_change, Q, I):
        """
        Initialize this QuotientMap.
        """
        self.__M_OK_change = M_OK_change
        self.__Q = Q
        self.__K = K
        self.__I = I
        self.__L, self.__from_L, self.__to_L = K.absolute_vector_space()

    def __call__(self, x):
        """
        Apply this QuotientMap to an element of the number field.

        INPUT:
            x -- an element of the field
        """
        v = self.__to_L(x)
        w = v * self.__M_OK_change
        return self.__Q( list(w) )

    def __repr__(self):
        """
        Return a string representation of this QuotientMap.
        """
        return "Partially defined quotient map from %s to an explicit vector space representation for the quotient of the ring of integers by (p,I) for the ideal I=%s."%(self.__K, self.__I)

class LiftMap:
    """
    Class to hold data needed by lifting maps from residue fields to
    number field orders.
    """
    def __init__(self, OK, M_OK_map, Q, I):
        """
        Initialize this LiftMap.
        """
        self.__I = I
        self.__OK = OK
        self.__Q = Q
        self.__M_OK_map = M_OK_map

    def __call__(self, x):
        """
        Apply this LiftMap to an element of the residue field.
        """
        # This lifts to OK tensor F_p
        v = self.__Q.lift(x)
        # This lifts to ZZ^n (= OK)
        w = v.lift()
        # Write back in terms of K
        z = w * self.__M_OK_map
        return self.__OK(z.list())

    def __repr__(self):
        """
        Return a string representation of this QuotientMap.
        """
        return "Lifting map to %s from quotient of integers by %s"%(self.__OK, self.__I)

def quotient_char_p(I, p):
    """
    Given an integral ideal I that contains a prime number p, compute
    a vector space V = (OK mod p) / (I mod p), along with a
    homomorphism OK --> V and a section V --> OK.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_ideal import quotient_char_p

        sage: K.<i> = NumberField(x^2 + 1); O = K.maximal_order(); I = K.fractional_ideal(15)
        sage: quotient_char_p(I, 5)[0]
        Vector space quotient V/W of dimension 2 over Finite Field of size 5 where
        V: Vector space of dimension 2 over Finite Field of size 5
        W: Vector space of degree 2 and dimension 0 over Finite Field of size 5
        Basis matrix:
        []
        sage: quotient_char_p(I, 3)[0]
        Vector space quotient V/W of dimension 2 over Finite Field of size 3 where
        V: Vector space of dimension 2 over Finite Field of size 3
        W: Vector space of degree 2 and dimension 0 over Finite Field of size 3
        Basis matrix:
        []

        sage: I = K.factor(13)[0][0]; I
        Fractional ideal (-3*i - 2)
        sage: I.residue_class_degree()
        1
        sage: quotient_char_p(I, 13)[0]
        Vector space quotient V/W of dimension 1 over Finite Field of size 13 where
        V: Vector space of dimension 2 over Finite Field of size 13
        W: Vector space of degree 2 and dimension 1 over Finite Field of size 13
        Basis matrix:
        [1 8]
    """
    if not I.is_integral():
        raise ValueError, "I must be an integral ideal."

    K    = I.number_field()
    OK   = K.maximal_order()  # will in the long run only really need a p-maximal order.
    M_OK = OK.free_module()
    M_I  = I.free_module()

    # Now we have to quite explicitly find a way to compute
    # with OK / I viewed as a quotient of two F_p vector spaces,
    # and itself viewed as an F_p vector space.

    # Step 1. Write each basis vector for I (as a ZZ-module)
    # in terms of the basis for OK.

    B_I = M_I.basis()
    M_OK_mat = M_OK.basis_matrix()
    M_OK_change = M_OK_mat**(-1)
    B_I_in_terms_of_M = M_I.basis_matrix() * M_OK_change

    # Step 2. Define "M_OK mod p" to just be (F_p)^n and
    # "M_I mod p" to be the reduction mod p of the elements
    # compute in step 1.

    n = K.degree()
    k = FiniteField(p)
    M_OK_modp = k**n
    B_mod = B_I_in_terms_of_M.change_ring(k)
    M_I_modp = M_OK_modp.span(B_mod.row_space())

    # Step 3. Compute the quotient of these two F_p vector space.

    Q = M_OK_modp.quotient(M_I_modp)

    # Step 4. Now we get the maps we need from the above data.

    K_to_Q = QuotientMap(K, M_OK_change, Q, I)
    Q_to_OK = LiftMap(OK, M_OK_mat, Q, I)

    return Q, K_to_Q, Q_to_OK



