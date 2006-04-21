"""
Finite Fields


EXAMPLES:
Finite Fields support iteration, starting with 0.
    sage: k = GF(9, 'a')
    sage: i = 0
    sage: for x in k: print i, x; i+=1
    0 0
    1 1
    2 2
    3 a
    4 a + 1
    5 a + 2
    6 2*a
    7 2*a + 1
    8 2*a + 2
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

import random, weakref

import arith
import field
import polynomial_ring
import sage.databases.conway
import sage.misc.defaults
import polynomial_element
import multi_polynomial_element

import integer
import rational

import sage.libs.pari.all as pari

import finite_field_element
import integer_mod_ring

from sage.structure.element import RingElement

import sage.interfaces.all

_objsFiniteField = {}
def FiniteField(order, name='a'):
    """
    Return a finite field of given order with generator labeled by the given name.

    INPUT:
        order -- int
        name -- string (default: 'a')

    EXAMPLES:
        sage: k, a = GF(9).objgen()
        sage: k
        Finite Field in a of size 3^2
        sage: k.assign_names(['b'])
        sage: k
        Finite Field in b of size 3^2
        sage: GF(9,'a')
        Finite Field in a of size 3^2
    """
    global _objsFiniteField
    key = (order, name)
    if _objsFiniteField.has_key(key):
        x = _objsFiniteField[key]()
        if x != None and x.variable_name() == name:      # see example above for why this is necessary
            return x
    if arith.is_prime(order):
        R = FiniteField_prime_modn(order, name)
    else:
        R = FiniteField_ext_pari(order, name)
    _objsFiniteField[key] = weakref.ref(R)
    return R

def GF(order, name='a'):
    """
    Synonym for FiniteField.
    """
    return FiniteField(order, name)

def is_FiniteField(x):
    return isinstance(x, FiniteField_generic)

def is_PrimeFiniteField(x):
    return isinstance(x, FiniteField_prime_modn)


class FiniteField_generic(field.Field):
    def __init__(self):
        """
        EXAMPLES:
            sage: K = GF(7); K
            Finite Field of size 7
            sage: loads(K.dumps()) == K
            True
            sage: GF(7^10)
            Finite Field in a of size 7^10
            sage: K = GF(7^10, 'a'); K
            Finite Field in a of size 7^10
            sage: loads(K.dumps()) == K
            True
        """
        raise NotImplementedError

    def _latex_(self):
        if self.degree() > 1:
            e = "^{%s}"%self.degree()
        else:
            e = ""
        return "\\mbox{\\rm F}_{%s%s}"%(self.characteristic(), e)

    def _gap_init_(self):
        return 'GF(%s)'%self.order()

    def __cmp__(self, other):
        """
        Compares this finite field with other.  Two finite fields are
        equal if and only if they have the same cardinality *and* the
        defining polynomials are the same.

        EXAMPLES:
            sage: FiniteField(3**2) == FiniteField(3**3)
            False
            sage: FiniteField(3**2) == FiniteField(3**2)
            True
            sage: FiniteField(3**2,'beta') == FiniteField(3**2,'alpha')
            False
            sage: FiniteField(3**2,'beta') == FiniteField(3**2,'beta')
            True
        """
        if self is other: return 0
        if not isinstance(other, FiniteField_generic):
            return -1
        if self.characteristic() < other.characteristic():
            return -1
        elif self.characteristic() > other.characteristic():
            return 1
        if self.variable_name() != other.variable_name():
            return -1
        if self.order() < other.order():
            return -1
        elif self.order()== other.order() and \
                 (self.degree() == 1 or self.polynomial() == other.polynomial()):
            return 0
        return 1

##     def __getitem__(self, n):
##         """
##         Returns $n$-th element of the field.  The ordering is
##         not randomized (though it could conceivably change from
##         one version of SAGE to another).

##         EXAMPLES:
##             sage: k = GF(8, 'a')
##             sage: k[0]
##             0
##             sage: k[1]
##             1
##             sage: k[7]
##             a^2 + a + 1
##         """
##         if n < 0 or n >= self.order():
##             raise IndexError, "n (=%s) must be between 0 and the order %s of the field."%(\
##                 n, self.order())
##         V = self.vector_space()
##         return self(V[n])

    def __iter__(self):
        for v in self.vector_space():
            yield self(v)

    def gen(self):
        raise NotImplementedError

    def zeta_order(self):
        return self.multiplicative_generator().multiplicative_order()

    def zeta(self, n=None):
        """
        Returns an element of multiplicative order n in this this
        finite field, if there is one.  Raises a ValueError if there
        is not.

        EXAMPLES:
            sage: k = GF(7)
            sage: k.zeta()
            3
            sage: k.zeta().multiplicative_order()
            6
            sage: k.zeta(3)
            2
            sage: k.zeta(3).multiplicative_order()
            3
            sage: k = GF(49)
            sage: k.zeta().multiplicative_order()
            48
            sage: k.zeta(6)
            3
        """
        z = self.multiplicative_generator()
        if n is None:
            return z
        else:
            n = integer.Integer(n)
            m = z.multiplicative_order()
            if m % n != 0:
                raise ValueError, "No %sth root of unity in self"%n
            return z**(m//n)

    def multiplicative_generator(self):
        """
        Return a generator for the multiplicative group of this field.
        The generator is not randomized, though it could change from
        one version of SAGE to another.

        EXAMPLES:
            sage: k = GF(997)
            sage: k.multiplicative_generator()
            7
            sage: k = GF(11**3, name='a')
            sage: k.multiplicative_generator()
            a
        """
        try:
            return self.__multiplicative_generator
        except AttributeError:
            if self.degree() == 1:
                self.__multiplicative_generator = self(arith.primitive_root(self.order()))
                return self.__multiplicative_generator
            n = self.order() - 1
            a = self.gen(0)
            if a.multiplicative_order() == n:
                self.__multiplicative_generator = a
                return a
            for a in self:
                if a == 0:
                    continue
                if a.multiplicative_order() == n:
                    self.__multiplicative_generator = a
                    return a

    def ngens(self):
        """
        The number of generators of the finite field.  Always 1.

        EXAMPLES:
            sage: k = FiniteField(3**4)
            sage: k.ngens()
            1
        """
        return 1

    def is_field(self):
        """
        Returns whether or not the finite field is a field, i.e.,
        always returns True.

        EXAMPLES:
            sage: k = FiniteField(3**4)
            sage: k.is_field()
            True
        """
        return True

    def is_finite(self):
        return True

    def order(self):
        raise NotImplementedError

    def cardinality(self):
        """
        Same as self.order().
        """
        return self.order()

    def unit_group_exponent(self):
        """
        The exponent of the unit group of the finite field.  For a
        finite field, this is always the order minus 1.

        EXAMPLES:
            sage: k = GF(2**10)
            sage: k.order()
            1024
            sage: k.unit_group_exponent()
            1023
        """
        return self.order() - 1


    def random_element(self, bound=None):
        """
        A random element of the finite field.

        INPUT:
            bound -- ignored

        EXAMPLES:
            sage.: k = GF(2**10, 'a')
            sage.: k.random_element()
            a^9 + a
        """
        if self.degree() == 1:
            return self(random.randrange(self.order()))
        v = self.vector_space().random_element()
        return self(v)

    def polynomial(self):
        raise NotImplementedError

    def polynomial_ring(self):
        """
        Returns the polynomial ring over the prime subfield in the
        same variable as this finite field.

        EXAMPLES:
            sage: k = FiniteField(3**4, "alpha")
            sage: k.polynomial_ring()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
        """
        try:
            return self.__polynomial_ring
        except AttributeError:
            self.__polynomial_ring = polynomial_ring.PolynomialRing(
                FiniteField(self.characteristic()), self.variable_name())
            return self.__polynomial_ring

    def vector_space(self):
        try:
            return self.__vector_space
        except AttributeError:
            import sage.modules.all
            V = sage.modules.all.VectorSpace(self.prime_subfield(),self.degree())
            self.__vector_space = V
            return V

class FiniteField_ext_pari(FiniteField_generic):
    """
    Finite Field of order q, where q is a nontrivial prime power.
    (Implemented using PARI mod's.)

    Create with the command
          FiniteField(order)

    INPUT:
        q -- int, prime power, order of the finite field
        name -- string (default: 'a'), string for printing the generator

    OUTPUT:
        FiniteField -- finite field of order q.

    EXAMPLES:
        sage: k = FiniteField(9, 'a')
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
        sage: a.charpoly()
        x^2 + 2*x + 2
        sage: [a**i for i in range(8)]
        [1, a, a + 1, 2*a + 1, 2, 2*a, 2*a + 2, a + 2]

    Fields can be coerced into sets or list and iterated over:
        sage: list(k)
        [0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2]

    The following is a native Python set:
        sage: set(k)
        set([a, 2*a + 1, a + 2, a + 1, 2*a + 2, 1, 0, 2, 2*a])

    And the following is a SAGE enumerated set:
        sage: EnumeratedSet(k)
        {a, 2*a + 1, a + 2, a + 1, 2*a + 2, 1, 0, 2, 2*a}

    We can also make a list via comprehension:
        sage: [x for x in k]
        [0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2]

    Next we compute with the finite field of order 16, where
    the name is named b.
        sage: k16 = FiniteField(16, "b")
        sage: z = k16.gen()
        sage: z
        b
        sage: z.charpoly()
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
        sage: K = GF(7)
        sage: loads(K.dumps()) == K
        True
        sage: K = GF(7^10)
        sage: loads(K.dumps()) == K
        True
        sage: K = GF(7^10, 'a')
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

        sage: K = GF(10007^10)
        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, q, name='a', modulus=None):
        """
        Create finite field of order q with variable printed as name.

        INPUT:
            q -- integer, size of the finite field, not prime
            name -- (optional:default 'a') variable used for
                    printing element of the finite field.  Also,
                    two finite fields are considered equal
                    if they have the same variable name, and not otherwise.
        OUTPUT:
            FiniteField -- a finite field of order q with given variable name.

        EXAMPLES:
            sage: FiniteField(17)
            Finite Field of size 17
            sage: FiniteField(2^10)
            Finite Field in a of size 2^10
            sage: FiniteField(3^5, "b")
            Finite Field in b of size 3^5
            sage: FiniteField(3^5, "b").gen()
            b

        You can also create a finite field using GF, which is a synonym
        for FiniteField.
            sage: GF(19**2)
            Finite Field in a of size 19^2
        """
        q = integer.Integer(q)
        if q < 2:
            raise ArithmeticError, "q (=%s) must be a prime power"%q
        F = q.factor()
        if len(F) > 1:
            raise ArithmeticError, "q (=%s) must be a prime power"%q
        self.assign_names(name)
        self.__char = F[0][0]
        self.__pari_one = pari.pari(1).Mod(self.__char)
        self.__degree = F[0][1]
        if F[0][1] <= 1:
            raise ValueError, "The size of the finite field must not be prime."
        self.__order = q
        self.__is_field = True
        if modulus is None:
            if exists_conway_polynomial(self.__char, self.__degree):
                modulus = conway_polynomial(self.__char, self.__degree)
            else:
                # The following is fast/deterministic, but has serious problems since
                # it crashes on 64-bit machines, and I can't figure out why:
                #     self.__pari_modulus = pari.pari.finitefield_init(self.__char, self.__degree, self.variable_name())
                # So instead we iterate through random polys until we find an irreducible one.

                R = polynomial_ring.PolynomialRing(GF(self.__char))
                while True:
                    modulus = R.random_element(self.__degree)
                    modulus = modulus.monic()
                    if modulus.degree() == self.__degree and modulus.is_irreducible():
                        break
        assert not (modulus is None)
        self.__modulus = modulus
        f = pari.pari(str(modulus))
        self.__pari_modulus = f.subst('x', 'a') * self.__pari_one
        self.__gen = finite_field_element.FiniteFieldElement(self,
            pari.pari('a'))


    def __cmp__(self, other):
        if not isinstance(other, FiniteField_ext_pari):
            return -1
        if (self is other) or (self.__order == other.__order and
                               self.variable_name() == other.variable_name() \
                               and self.__modulus == other.__modulus):
            return 0
        return 1

    def _pari_one(self):
        """
        The PARI object Mod(1,p).  This is implementation specific
        and should be ignored by users.
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
            sage: GF(19**2, 'a')._pari_modulus()
            Mod(1, 19)*a^2 + Mod(18, 19)*a + Mod(2, 19)

            sage: GF(13**3, 'a')._pari_modulus()
            Mod(1, 13)*a^3 + Mod(2, 13)*a + Mod(11, 13)

        Note that the PARI modulus is always in terms of a, even if
        the field variable isn't.  This is because the specific choice
        of variable name has meaning in PARI, i.e., it can't be
        arbitrary.
            sage: FiniteField(2**4, "b")._pari_modulus()
            Mod(1, 2)*a^4 + Mod(1, 2)*a + Mod(1, 2)
        """
        return self.__pari_modulus

    def is_prime_field(self):
        return False

    def is_prime(self):
        return False

    def gen(self, n=0):
        """
        Return chosen generator of the finite field.  This generator
        is a root of the defining polynomial of the finite field, and
        is guaranteed to be a generator for the multiplicative group.

        INPUT:
            nothing

        OUTPUT:
            FiniteFieldElement -- field generator of finite field

        EXAMPLES:
            sage: FiniteField(2**4, "b").gen()
            b
            sage: k = FiniteField(3**4, "alpha")
            sage: a = k.gen()
            sage: a
            alpha
            sage: a**4
            alpha^3 + 1
        """
        return self.__gen

    def characteristic(self):
        """
        Returns the characteristic of the finite field, which is a
        prime int.

        EXAMPLES:
            sage: k = FiniteField(3**4)
            sage: k.characteristic()
            3
        """
        return self.__char

    def modulus(self):
        return self.__modulus

    def degree(self):
        """
        Returns the degree of the finite field, which is a positive
        integer.

        EXAMPLES:
            sage: FiniteField(3).degree()
            1
            sage: FiniteField(3**20).degree()
            20
        """
        return self.__degree

    def __repr__(self):
        return "Finite Field in %s of size %s^%s"%(self.variable_name(), self.__char, self.__degree)

    def __call__(self, x):
        r"""
        Coerce x into the finite field.

        INPUT:
            x -- object

        OUTPUT:
            FiniteFieldElement -- if possible, makes a finite field element from x.

        EXAMPLES:
            sage: k = GF(3^4)
            sage: b = k(5)
            sage: b.parent()
            Finite Field in a of size 3^4
            sage: a = k.gen()
            sage: k(a + 2)
            a + 2

        Constant polynomials coerce into finite fields:
            sage: R = QQ['x']
            sage: k, a = GF(5^2).objgen()
            sage: k(R(2/3))
            4
            sage: R, x = k['x'].objgen()
            sage: k(R(3))
            3

        Nonconstant polynomials do not coerce:
            sage: k(x)
            Traceback (most recent call last):
            ...
            TypeError: no coercion of non-constant polynomial x into Finite Field in a of size 5^2 defined.
            sage: k(R(a))
            a

        Multivariate polynomials also coerce:
            sage: R = k['x,y,z']; R
            Polynomial Ring in x, y, z over Finite Field in a of size 5^2
            sage: k(R(2))
            2
            sage: R = QQ['x,y,z']
            sage: k(R(1/5))
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce 1/5 into Finite Field in a of size 5^2.


        \note{Finite Fields are currently implemented using
        polynomials modulo p and the PARI ffinit function.  In
        particular, we do not use Conway polynomials and do \emph{not}
        yet define natural consistent inclusion maps between different
        finite fields.}

        Gap elements can also be coerced into finite fields.

            sage: F = GF(8, 'a')
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

            sage: F = GF(16)
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
        FiniteFieldElements.

        AUTHOR:
            -- David Joyner (2005-11)
            -- Martin Albrecht (2006-01-23)
            -- Martin Albrecth (2006-03-06): added coercion from string
        """
        if isinstance(x, finite_field_element.FiniteFieldElement):
            if x.parent() == self:
                return x
            else:
                # This is where we *would* do coercion from one finite field to another...
                raise TypeError, "No coercion of %s into %s defined."%(x, self)

        elif sage.interfaces.all.is_GapElement(x):
            try:
                return gap_to_sage(x, self)
            except (ValueError, IndexError, TypeError):
                raise TypeError, "error coercing %s into %s"%(x, self)

        if isinstance(x, (int, long, integer.Integer, rational.Rational,
                          pari.pari_gen)):

            return finite_field_element.FiniteFieldElement(self, x)

        elif isinstance(x, (multi_polynomial_element.MPolynomial, polynomial_element.Polynomial)):
            if x.is_constant():
                return self(x.constant_coefficient())
            else:
                raise TypeError, "no coercion of non-constant polynomial %s into %s defined."%(x,self)

        elif isinstance(x, str):
            x = pari.pari(x)
            t = x.type()
            if t == 't_POL':
                if (x.variable() == self.variable_name() \
                    and x.polcoeff(0).type()[2] == 'I'): #t_INT and t_INTMOD
                    return self(x)
            if t[2] == 'I': #t_INT and t_INTMOD
                return self(x)
            raise TypeError, "string element does not match this finite field"

        try:
            if x.parent() == self.vector_space():
                name = self.variable_name()
                x = pari.pari('+'.join(['%s*%s^%s'%(x[i],name, i) for i in range(self.degree())]))
                return finite_field_element.FiniteFieldElement(self, x)
        except AttributeError:
            pass
        try:
            return finite_field_element.FiniteFieldElement(self, integer.Integer(x))
        except TypeError:
            raise TypeError, "no coercion of %s into %s defined."%(x, self)

    def _coerce_(self, x):
        if isinstance(x, (int, long, integer.Integer)):
            return self(x)
        if isinstance(x, RingElement):
            K = x.parent()
            if K is self:
                return x
            if K.characteristic() == self.characteristic():
                if K.degree() == 1:
                    return self(int(x))
                elif self.degree() % K.degree() == 0:
                    # This is where we *would* do coercion from one nontrivial finite field to another...
                    raise NotImplementedError, "nontrivial finite field coercions not implemented"
        raise TypeError

    def __len__(self):
        """
        The number of elements of the finite field.

        EXAMPLES:
            sage: k = GF(2**10)
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
            sage: k = GF(2**10)
            sage: k
            Finite Field in a of size 2^10
            sage: k.order()
            1024
        """
        return self.__order

    def polynomial(self):
        """
        Return the irreducible characteristic polynomial of the
        generator of this finite field, i.e., the polynomial f(x) so
        elements of the finite field as elements modulo f.

        EXAMPLES:
            sage: k = FiniteField(17)
            sage: k.polynomial()
            x
            sage: k = FiniteField(9)
            sage: k.polynomial()
            x^2 + 2*x + 2
        """
        try:
            return self.__polynomial
        except  AttributeError:
            R = polynomial_ring.PolynomialRing(FiniteField(self.characteristic()))
            self.__polynomial = R(self._pari_modulus())
        return self.__polynomial


class FiniteField_prime_modn(FiniteField_generic, integer_mod_ring.IntegerModRing_field):
    def __init__(self, p, name=None):
        p = integer.Integer(p)
        if not arith.is_prime(p):
            raise ArithmeticError, "p (=%s) must be prime."%p
        integer_mod_ring.IntegerModRing_field.__init__(self, p)
        self.__char = p
        self.__gen = self(1)  # self(int(pari.pari(p).znprimroot().lift()))
        self.assign_names(name)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def _coerce_(self, x):
        if isinstance(x, (int, long, integer.Integer)):
            return self(x)
        if isinstance(x, RingElement):
            K = x.parent()
            if K is self:
                return x
        raise TypeError, "Unable to coerce %s into %s"%(x, self)

    def characteristic(self):
        return self.__char

    def modulus(self):
        try:
            return self.__modulus
        except AttributeError:
            x = polynomial_ring.PolynomialRing(self).gen()
            self.__modulus = x - 1
        return self.__modulus

    def is_prime_field(self):
        return True

    def is_prime(self):
        return True

    def polynomial(self):
        try:
            return self.__polynomial
        except  AttributeError:
            self.__polynomial = polynomial_ring.PolynomialRing(self)([0,1])
            return self.__polynomial

    def order(self):
        return self.__char

    def gen(self, n=0):
        return self.__gen

    def __repr__(self):
        return "Finite Field of size %s"%(self.order())

    def __iter__(self):
        for i in xrange(self.order()):
            yield self(i)

    def degree(self):
        """
        Returns the degree of the finite field, which is a positive
        integer.

        EXAMPLES:
            sage: FiniteField(3).degree()
            1
            sage: FiniteField(3**20).degree()
            20
        """
        return 1

##     def __getitem__(self, n):
##         if n < 0 or n >= self.modulus():
##             raise IndexError
##         return self(n)


##################################################################

def conway_polynomial(p, n):
    """
    Return the Conway polynomial of degree n over GF(p), which is
    loaded from a table.

    If the requested polynomial is not known, this function raises a
    RuntimeError exception.

    INPUT:
        p -- int
        n -- int

    OUTPUT:
        Polynomial -- a polynomial over the prime finite field GF(p).

    NOTE: The first time this function is called a table is read from
    disk, which takes a fraction of a second.  Subsequent calls do not
    require reloading the table.

    See also the ConwayPolynomials() object, which is a table of
    Conway polynomials.   For example, if c=ConwayPolynomials, then
    c.primes() is a list of all primes for which the polynomials are
    known, and for a given prime p,  c.degree(p) is a list of all
    degrees for which the Conway polynomials are known.

    EXAMPLES:
        sage: conway_polynomial(2,5)
        x^5 + x^2 + 1
        sage: conway_polynomial(101,5)
        x^5 + 2*x + 99
        sage: conway_polynomial(97,101)
        Traceback (most recent call last):
        ...
        RuntimeError: Conway polynomial over F_97 of degree 101 not in database.
    """
    (p,n)=(int(p),int(n))
    R = polynomial_ring.PolynomialRing(GF(p), 'x')
    try:
        return R(sage.databases.conway.ConwayPolynomials()[p][n])
    except KeyError:
        raise RuntimeError, \
              "Conway polynomial over F_%s of degree %s not in database."%(p,n)

def exists_conway_polynomial(p, n):
    """
    Return True if the Conway polynomial over F_p of degree n is in the
    database and False otherwise.

    If the Conway polynomial is in the database, to obtain it use the
    command conway_polynomial(p,n).

    EXAMPLES:
        sage: exists_conway_polynomial(2,3)
        True
        sage: exists_conway_polynomial(2,-1)
        False
        sage: exists_conway_polynomial(97,200)
        False
        sage: exists_conway_polynomial(6,6)
        False
    """
    return sage.databases.conway.ConwayPolynomials().has_polynomial(p,n)


def gap_to_sage(x, F):
    """
    INPUT:
        x -- gap finite field element
        F -- SAGE finite field
    OUTPUT:
        element of F

    EXAMPLES:
        sage: x = gap('Z(13)')
        sage: F = GF(13)
        sage: F(x)
        2
        sage: F(gap('0*Z(13)'))
        0
        sage: F = GF(13^2)
        sage: x = gap('Z(13)')
        sage: F(x)
        2
        sage: x = gap('Z(13^2)^3')
        sage: F(x)
        12*a + 11
        sage: F.multiplicative_generator()^3
        12*a + 11

    AUTHOR:
        -- David Joyner and William Stein
    """
    s = str(x)
    if s[:2] == '0*':
        return F(0)
    i1 = s.index("(")
    i2 = s.index(")")
    q  = eval(s[i1+1:i2].replace('^','**'))
    if q == F.order():
        K = F
    else:
        K = FiniteField(q)
    if s.find(')^') == -1:
        e = 1
    else:
        e = int(s[i2+2:])
    if F.degree() == 1:
        g = int(sage.interfaces.all.gap.eval('Int(Z(%s))'%q))
    else:
        g = K.multiplicative_generator()
    return F(K(g**e))

