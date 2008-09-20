"""
Finite Extension Fields implemented via PARI.
"""

#*****************************************************************************
#       Copyright (C) 2005,2007 William Stein <wstein@gmail.com>
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

import polynomial.polynomial_element as polynomial_element
import polynomial.multi_polynomial_element as multi_polynomial_element

import integer
import rational
import integer_mod
import integer_mod_ring

import sage.libs.pari.all as pari

import finite_field_element

from sage.structure.element import RingElement
from sage.rings.ring import FiniteField as FiniteField_generic

from sage.structure.parent_gens import normalize_names, ParentWithGens

import sage.interfaces.gap

class FiniteField_ext_pari(FiniteField_generic):
    r"""
    Finite Field of order q, where q is a nontrivial prime power.
    (Implemented using PARI mod's.) This implementation is the default
    implementation for $q \geq 2^{16}$.

    Create with the command
          FiniteField(order)

    INPUT:
        q -- int, prime power, order of the finite field
        name -- string (default: 'a'), string for printing the generator

    OUTPUT:
        FiniteField_ext_pari -- finite field of order q.

    EXAMPLES:
        sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
        sage: k = FiniteField_ext_pari(9, 'a')
        sage: k
        Finite Field in a of size 3^2
        sage: k.is_field()
        True
        sage: k.characteristic()
        3
        sage: a = k.gen()
        sage: a
        a
        sage: a.parent()
        Finite Field in a of size 3^2
        sage: a.charpoly('x')
        x^2 + 2*x + 2
        sage: [a^i for i in range(8)]
        [1, a, a + 1, 2*a + 1, 2, 2*a, 2*a + 2, a + 2]

    Fields can be coerced into sets or list and iterated over:
        sage: list(k)
        [0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2]

    The following is a native Python set:
        sage: set(k)
        set([0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2])

    And the following is a SAGE enumerated set:
        sage: EnumeratedSet(k)
        {0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2}

        We can also make a list via comprehension:
        sage: [x for x in k]
        [0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2]

    Next we compute with the finite field of order 16, where
    the name is named b.
        sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
        sage: k16 = FiniteField_ext_pari(16, "b")
        sage: z = k16.gen()
        sage: z
        b
        sage: z.charpoly('x')
        x^4 + x + 1
        sage: k16.is_field()
        True
        sage: k16.characteristic()
        2
        sage: z.multiplicative_order()
        15

    Of course one can also make prime finite fields.
        sage: k = FiniteField(7)

    Note that the generator is 1:
        sage: k.gen()
        1
        sage: k.gen().multiplicative_order()
        1

    Illustration of dumping and loading:
        sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
        sage: K = FiniteField(7)
        sage: loads(K.dumps()) == K
        True
        sage: K = FiniteField_ext_pari(7^10, 'b')
        sage: loads(K.dumps()) == K
        True
        sage: K = FiniteField_ext_pari(7^10, 'a')
        sage: loads(K.dumps()) == K
        True

    In this example $K$ is large enough that Conway polynomials
    are not used.  Note that when the field is dumped the defining
    polynomial $f$ is also dumped.  Since $f$ is determined by a
    random algorithm, it's important that $f$ is dumped as part of
    $K$.  If you quit \sage and restart and remake a finite field
    of the same order (and the order is large enough so that there
    is no Conway polynomial), then defining polynomial is probably
    different.  However, if you load a previously saved field, that
    will have the same defining polynomial.

        sage: K = GF(10007^10, 'a')
        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, q, name, modulus=None):
        """
        Create finite field of order q with variable printed as name.

        INPUT:
            q -- integer, size of the finite field, not prime
            name -- variable used for printing element of the finite
                    field.  Also, two finite fields are considered
                    equal if they have the same variable name, and not
                    otherwise.
            modulus -- you may provide a minimal polynomial to use for
                       reduction or None to force a random or conway
                       irreducible polynomial. (default: None, a conway
                       polynomial is used if found. Otherwise a random
                       polynomial is used)

        OUTPUT:
            FiniteField_ext_pari -- a finite field of order q with given variable name.

        EXAMPLES:
            sage: FiniteField(17)
            Finite Field of size 17
            sage: FiniteField(2^10, 'c')
            Finite Field in c of size 2^10
            sage: FiniteField(3^5, "b")
            Finite Field in b of size 3^5
            sage: FiniteField(3^5, "b").gen()
            b

        You can also create a finite field using GF, which is a synonym
        for FiniteField.
            sage: GF(19^2, 'a')
            Finite Field in a of size 19^2
        """
        from finite_field import FiniteField as GF
        q = integer.Integer(q)
        if q < 2:
            raise ArithmeticError, "q must be a prime power"
        F = q.factor()
        if len(F) != 1:
            raise ArithmeticError, "q must be a prime power"

        if F[0][1] > 1:
            base_ring = GF(F[0][0])
        else:
            raise ValueError, "The size of the finite field must not be prime."
            #base_ring = self

        ParentWithGens.__init__(self, base_ring, name, normalize=True)

        self._kwargs = {}
        self.__char = F[0][0]
        self.__pari_one = pari.pari(1).Mod(self.__char)
        self.__degree = F[0][1]
        self.__order = q
        self.__is_field = True
        if modulus is None:
            from finite_field import conway_polynomial
            from finite_field import exists_conway_polynomial

            if exists_conway_polynomial(self.__char, self.__degree):
                modulus = conway_polynomial(self.__char, self.__degree)
            else:
                # The following is fast/deterministic, but has serious problems since
                # it crashes on 64-bit machines, and I can't figure out why:
                #     self.__pari_modulus = pari.pari.finitefield_init(self.__char, self.__degree, self.variable_name())
                # So instead we iterate through random polys until we find an irreducible one.

                R = GF(self.__char)['x']
                while True:
                    modulus = R.random_element(self.__degree)
                    modulus = modulus.monic()
                    if modulus.degree() == self.__degree and modulus.is_irreducible():
                        break
        assert not (modulus is None)
        if isinstance(modulus, (list, tuple)):
            modulus = GF(self.__char)['x'](modulus)
        self.__modulus = modulus
        f = pari.pari(str(modulus))
        self.__pari_modulus = f.subst(modulus.parent().variable_name(), 'a') * self.__pari_one
        self.__gen = finite_field_element.FiniteField_ext_pariElement(self, pari.pari('a'))

        self._zero_element = self(0)
        self._one_element = self(1)

    def __cmp__(self, other):
        """
        EXAMPLE:
            sage: k = GF(7^20,'a')
            sage: k == loads(dumps(k))
            True
        """
        if not isinstance(other, FiniteField_ext_pari):
            return cmp(type(self), type(other))
        return cmp((self.__order, self.variable_name()), (other.__order, other.variable_name()))

    def _pari_one(self):
        r"""
        The \PARI object Mod(1,p).  This is implementation specific
        and should be ignored by users.

        EXAMPLE:
            sage: k = GF(7^20,'a')
            sage: k._pari_one()
            Mod(1, 7)
        """
        return self.__pari_one

    def _pari_modulus(self):
        """
        The polynomial mod p that defines the finite field, as a PARI
        object.  This is implementation specific, and some finite fields
        might not be implemented using PARI, so you should avoid using
        this function.

        INPUT:  nothing
        OUTPUT:
            gen -- a pari polynomial gen
        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(19^2, 'a')._pari_modulus()
            Mod(1, 19)*a^2 + Mod(18, 19)*a + Mod(2, 19)

            sage: FiniteField_ext_pari(13^3, 'a')._pari_modulus()
            Mod(1, 13)*a^3 + Mod(2, 13)*a + Mod(11, 13)

        Note that the PARI modulus is always in terms of a, even if
        the field variable isn't.  This is because the specific choice
        of variable name has meaning in PARI, i.e., it can't be
        arbitrary.
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(2^4, "b")._pari_modulus()
            Mod(1, 2)*a^4 + Mod(1, 2)*a + Mod(1, 2)
        """
        return self.__pari_modulus

    def gen(self, n=0):
        """
        Return chosen generator of the finite field.  This generator
        is a root of the defining polynomial of the finite field.

        WARNING: The generator is not guaranteed to be a generator for
            the multiplicative group.  To obtain the latter, use
            multiplicative_generator().  Both gen() and
            multiplicative_generator() are random: the elements
            returned will in general differ between runs.

        INPUT:
            nothing

        OUTPUT:
            FiniteField_ext_pariElement -- field generator of finite field

        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(2^4, "b").gen()
            b
            sage: k = FiniteField_ext_pari(3^4, "alpha")
            sage: a = k.gen()
            sage: a
            alpha
            sage: a^4
            alpha^3 + 1

        """
        return self.__gen

    def characteristic(self):
        """
        Returns the characteristic of the finite field, which is a
        prime int.

        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(3^4, 'a')
            sage: k.characteristic()
            3
        """
        return self.__char

    def modulus(self):
        r"""
        Return the minimal polynomial of the generator of self in
        \code{self.polynomial_ring('x')}.

        EXAMPLES:
            sage: F.<a> = GF(7^20, 'a')
            sage: f = F.modulus(); f
            x^20 + x^12 + 6*x^11 + 2*x^10 + 5*x^9 + 2*x^8 + 3*x^7 + x^6 + 3*x^5 + 3*x^3 + x + 3

            sage: f(a)
            0
        """
        return self.__modulus

    def degree(self):
        """
        Returns the degree of the finite field, which is a positive
        integer.

        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField(3).degree()
            1
            sage: FiniteField_ext_pari(3^20, 'a').degree()
            20
        """
        return self.__degree

    def __call__(self, x):
        r"""
        Coerce x into the finite field.

        INPUT:
            x -- object

        OUTPUT:
            FiniteField_ext_pariElement -- if possible, makes a finite field element from x.

        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(3^4, 'a')
            sage: b = k(5)
            sage: b.parent()
            Finite Field in a of size 3^4
            sage: a = k.gen()
            sage: k(a + 2)
            a + 2

        Univariate polynomials coerce into finite fields by evaluating
        the polynomial at the field's generator:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: R.<x> = QQ[]
            sage: k, a = FiniteField_ext_pari(5^2, 'a').objgen()
            sage: k(R(2/3))
            4
            sage: k(x^2)
            a + 3
            sage: R.<x> = GF(5)[]
            sage: k(x^3-2*x+1)
            2*a + 4

            sage: x = polygen(QQ)
            sage: k(x^25)
            a

            sage: Q, q = FiniteField_ext_pari(5^7, 'q').objgen()
            sage: L = GF(5)
            sage: LL.<xx> = L[]
            sage: Q(xx^2 + 2*xx + 4)
            q^2 + 2*q + 4


        Multivariate polynomials only coerce if constant:
            sage: R = k['x,y,z']; R
            Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 5^2
            sage: k(R(2))
            2
            sage: R = QQ['x,y,z']
            sage: k(R(1/5))
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce


        \note{Finite Fields are currently implemented using
        polynomials modulo p and the PARI ffinit function.  In
        particular, we do not use Conway polynomials and do \emph{not}
        yet define natural consistent inclusion maps between different
        finite fields.}

        Gap elements can also be coerced into finite fields.

            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: F = FiniteField_ext_pari(8, 'a')
            sage: a = F.multiplicative_generator(); a
            a
            sage: b = gap(a^3); b
            Z(2^3)^3
            sage: F(b)
            a + 1
            sage: a^3
            a + 1

            sage: a = GF(13)(gap('0*Z(13)')); a
            0
            sage: a.parent()
            Finite Field of size 13

            sage: F = GF(16, 'a')
            sage: F(gap('Z(16)^3'))
            a^3
            sage: F(gap('Z(16)^2'))
            a^2

        You can also call a finite extension field with a string
        to produce an element of that field, like this:

            sage: k = GF(2^8, 'a')
            sage: k('a^200')
            a^4 + a^3 + a^2

        This is especially useful for fast conversions from Singular etc. to
        FiniteField_ext_pariElements.

        AUTHOR:
            -- David Joyner (2005-11)
            -- Martin Albrecht (2006-01-23)
            -- Martin Albrecth (2006-03-06): added coercion from string
        """
        if isinstance(x, finite_field_element.FiniteField_ext_pariElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                # canonically isomorphic finite fields
                return finite_field_element.FiniteField_ext_pariElement(self, x)
            else:
                # This is where we *would* do coercion from one finite field to another...
                raise TypeError, "no coercion defined"

        elif sage.interfaces.gap.is_GapElement(x):
            from sage.interfaces.gap import gfq_gap_to_sage
            try:
                return gfq_gap_to_sage(x, self)
            except (ValueError, IndexError, TypeError):
                raise TypeError, "no coercion defined"

        if isinstance(x, (int, long, integer.Integer, rational.Rational,
                          pari.pari_gen)):

            return finite_field_element.FiniteField_ext_pariElement(self, x)

        elif isinstance(x, multi_polynomial_element.MPolynomial):
            if x.is_constant():
                return self(x.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

        elif isinstance(x, polynomial_element.Polynomial):
            if x.is_constant():
                return self(x.constant_coefficient())
            else:
                return x.change_ring(self)(self.gen())

        elif isinstance(x, str):
            x = x.replace(self.variable_name(),'a')
            x = pari.pari(x)
            t = x.type()
            if t == 't_POL':
                if (x.variable() == 'a' \
                    and x.polcoeff(0).type()[2] == 'I'): #t_INT and t_INTMOD
                    return self(x)
            if t[2] == 'I': #t_INT and t_INTMOD
                return self(x)
            raise TypeError, "string element does not match this finite field"

        try:
            if x.parent() == self.vector_space():
                x = pari.pari('+'.join(['%s*a^%s'%(x[i], i) for i in range(self.degree())]))
                return finite_field_element.FiniteField_ext_pariElement(self, x)
        except AttributeError:
            pass
        try:
            return finite_field_element.FiniteField_ext_pariElement(self, integer.Integer(x))
        except TypeError, msg:
            raise TypeError, "%s\nno coercion defined"%msg

    def _coerce_impl(self, x):
        r"""
        Canonical coercion to \code{self}.

        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(4,'a')._coerce_(GF(2)(1)) # indirect doctest
            1
            sage: k = FiniteField_ext_pari(4,'a')
            sage: k._coerce_(k.0)
            a
            sage: FiniteField_ext_pari(4,'a')._coerce_(3)
            1
            sage: FiniteField_ext_pari(4,'a')._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion defined
            sage: FiniteField_ext_pari(8,'a')._coerce_(FiniteField_ext_pari(4,'a').0)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion defined
            sage: FiniteField_ext_pari(16,'a')._coerce_(FiniteField_ext_pari(4,'a').0)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion defined
            sage: k = FiniteField_ext_pari(8,'a')
            sage: k._coerce_(FiniteField(7,'a')(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion defined
        """
        if isinstance(x, (int, long, integer.Integer)):
            return self(x)

        if isinstance(x, finite_field_element.FiniteField_ext_pariElement) or integer_mod.is_IntegerMod(x):
            K = x.parent()
            if K is self:
                return x
            if isinstance(K, integer_mod_ring.IntegerModRing_generic) and K.characteristic() % self.characteristic() == 0:
                return self(int(x))
            if K.characteristic() == self.characteristic():
                if K.degree() == 1:
                    return self(int(x))
                elif self.degree() % K.degree() == 0:
                    # TODO: This is where we *would* do coercion from one nontrivial finite field to another...
                    raise TypeError, 'no canonical coercion defined'
        raise TypeError, 'no canonical coercion defined'

    def __len__(self):
        """
        The number of elements of the finite field.

        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(2^10, 'a')
            sage: k
            Finite Field in a of size 2^10
            sage: len(k)
            1024
        """
        return self.__order

    def order(self):
        """
        The number of elements of the finite field.

        EXAMPLES:
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(2^10,'a')
            sage: k
            Finite Field in a of size 2^10
            sage: k.order()
            1024
        """
        return self.__order

    def polynomial(self, name=None):
        """
        Return the irreducible characteristic polynomial of the
        generator of this finite field, i.e., the polynomial f(x) so
        elements of the finite field as elements modulo f.

        EXAMPLES:
            sage: k = FiniteField(17)
            sage: k.polynomial('x')
            x
            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(9,'a')
            sage: k.polynomial('x')
            x^2 + 2*x + 2
        """
        if name is None:
            name = self.variable_name()
        try:
            return self.__polynomial[name]
        except (AttributeError, KeyError):
            from finite_field import FiniteField as GF
            R = GF(self.characteristic())[name]
            f = R(self._pari_modulus())
            try:
                self.__polynomial[name] = f
            except (KeyError, AttributeError):
                self.__polynomial = {}
                self.__polynomial[name] = f
            return f

    def __hash__(self):
        """
        Return the hash of this field.

        EXAMPLES:
            sage: hash(GF(3,'b'))
            904200654                      # 32-bit
            -586939294780423730            # 64-bit
            sage: hash(GF(3,'a'))
            904200654                      # 32-bit
            -586939294780423730            # 64-bit
            sage: hash(GF(9,'a'))
            -417021630                     # 32-bit
            1006006598732398914            # 64-bit
            sage: hash(GF(9,'b'))
            995034271                      # 32-bit
            8600900932991911071            # 64-bit
        """
        return hash((self.__order, self.variable_name(), self.__modulus))

