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
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
import integer_mod

from sage.structure.element import RingElement
from sage.rings.ring import FiniteField as FiniteField_generic
from sage.rings.finite_field_givaro import FiniteField_givaro

import sage.interfaces.gap

from finite_field_c import FiniteField, is_FiniteField, is_PrimeFiniteField

import weakref

cache = {}

def FiniteField(order, name=None, modulus=None):
    """
    Return the globally unique finite field of given order with generator
    labeled by the given name and possibly with given modulus.

    INPUT:
        order --   int
        name --    string; must be specified in not a prime field
        modulus -- (optional) defining polynomial for field, i.e.,
                   generator of the field will be a root of this
                   polynomial; if not specified the choice of
                   definining polynomials can be arbitrary.

    EXAMPLES:
        sage: k = FiniteField(9, 'a'); k
        Finite Field in a of size 3^2
        sage: parent(a)
        Finite Field in a of size 3^2
        sage: charpoly(a, 'y')
        y^2 + 2*y + 2

    You can also use GF instead of FiniteField -- they are identical.
    """
    order = int(order)

    key = (order, name, modulus)
    if cache.has_key(key):
        return cache[key]
    # I have disabled weakref support for finite fields, because it isn't
    # really implemented in Pyrex.  - SEE track ticket #165
        #K = cache[key]()
        #if not K is None:
        #    return K

    if arith.is_prime(order):
        K = integer_mod_ring.IntegerModRing(order)
    else:
        if name is None:
            raise TypeError, "you must specify the generator name"
        if order < 2**16:   # todo -- re-enable
            K = FiniteField_givaro(order, name, modulus)
        else:
            K = FiniteField_ext_pari(order, name, modulus)

    #cache[key] = weakref.ref(K)
    #cache[key] = K
    return K


def is_FiniteField(x):
    return isinstance(x, FiniteField_generic)

def is_PrimeFiniteField(x):
    return isinstance(x, FiniteField_prime_modn)

GF = FiniteField

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
        sage: from sage.rings.finite_field import FiniteField_ext_pari
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
        sage: from sage.rings.finite_field import FiniteField_ext_pari
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
        sage: from sage.rings.finite_field import FiniteField_ext_pari
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
            sage: GF(19**2, 'a')
            Finite Field in a of size 19^2
        """
        q = integer.Integer(q)
        if q < 2:
            raise ArithmeticError, "q must be a prime power"
        F = q.factor()
        if len(F) > 1:
            raise ArithmeticError, "q must be a prime power"
        self._assign_names(name)
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

                R = polynomial_ring.PolynomialRing(GF(self.__char), 'x')
                while True:
                    modulus = R.random_element(self.__degree)
                    modulus = modulus.monic()
                    if modulus.degree() == self.__degree and modulus.is_irreducible():
                        break
        assert not (modulus is None)
        self.__modulus = modulus
        f = pari.pari(str(modulus))
        self.__pari_modulus = f.subst('x', 'a') * self.__pari_one
        self.__gen = finite_field_element.FiniteField_ext_pariElement(self, pari.pari('a'))


    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: GF(7)(2) == GF(7)(9)
            True
            sage: GF(7)(2) == GF(11)(2)
            False
            sage: GF(7)(2) == GF(8,'a')(2)
            False
            sage: GF(7)(2) == 2
            True
        """
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
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: FiniteField_ext_pari(19**2, 'a')._pari_modulus()
            Mod(1, 19)*a^2 + Mod(18, 19)*a + Mod(2, 19)

            sage: FiniteField_ext_pari(13**3, 'a')._pari_modulus()
            Mod(1, 13)*a^3 + Mod(2, 13)*a + Mod(11, 13)

        Note that the PARI modulus is always in terms of a, even if
        the field variable isn't.  This is because the specific choice
        of variable name has meaning in PARI, i.e., it can't be
        arbitrary.
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: FiniteField_ext_pari(2**4, "b")._pari_modulus()
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
            FiniteField_ext_pariElement -- field generator of finite field

        EXAMPLES:
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: FiniteField_ext_pari(2**4, "b").gen()
            b
            sage: k = FiniteField_ext_pari(3**4, "alpha")
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
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(3**4, 'a')
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
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: FiniteField(3).degree()
            1
            sage: FiniteField_ext_pari(3**20, 'a').degree()
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
            FiniteField_ext_pariElement -- if possible, makes a finite field element from x.

        EXAMPLES:
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(3^4, 'a')
            sage: b = k(5)
            sage: b.parent()
            Finite Field in a of size 3^4
            sage: a = k.gen()
            sage: k(a + 2)
            a + 2

        Constant polynomials coerce into finite fields:
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: R = QQ['x']
            sage: k, a = FiniteField_ext_pari(5^2, 'a').objgen()
            sage: k(R(2/3))
            4
            sage: R = k['x']
            sage: k(R(3))
            3

        Nonconstant polynomials do not coerce:
            sage: k(x)
            Traceback (most recent call last):
            ...
            TypeError: no coercion defined
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
            TypeError: unable to coerce


        \note{Finite Fields are currently implemented using
        polynomials modulo p and the PARI ffinit function.  In
        particular, we do not use Conway polynomials and do \emph{not}
        yet define natural consistent inclusion maps between different
        finite fields.}

        Gap elements can also be coerced into finite fields.

            sage: from sage.rings.finite_field import FiniteField_ext_pari
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
            if x.parent() == self:
                return x
            else:
                # This is where we *would* do coercion from one finite field to another...
                raise TypeError, "no coercion defined"

        elif sage.interfaces.gap.is_GapElement(x):
            try:
                return gap_to_sage(x, self)
            except (ValueError, IndexError, TypeError):
                raise TypeError, "no coercion defined"

        if isinstance(x, (int, long, integer.Integer, rational.Rational,
                          pari.pari_gen)):

            return finite_field_element.FiniteField_ext_pariElement(self, x)

        elif isinstance(x, (multi_polynomial_element.MPolynomial, polynomial_element.Polynomial)):
            if x.is_constant():
                return self(x.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

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

    def _coerce_(self, x):
        """
        Canonical coercion to self.

        EXAMPLES:
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: FiniteField_ext_pari(4,'a')._coerce_(GF(2)(1))
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

        if isinstance(x, (finite_field_element.FiniteField_ext_pariElement)) or integer_mod.is_IntegerMod(x):
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
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(2**10, 'a')
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
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(2**10,'a')
            sage: k
            Finite Field in a of size 2^10
            sage: k.order()
            1024
        """
        return self.__order

    def polynomial(self, name):
        """
        Return the irreducible characteristic polynomial of the
        generator of this finite field, i.e., the polynomial f(x) so
        elements of the finite field as elements modulo f.

        EXAMPLES:
            sage: k = FiniteField(17)
            sage: k.polynomial('x')
            x
            sage: from sage.rings.finite_field import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(9,'a')
            sage: k.polynomial('x')
            x^2 + 2*x + 2
        """
        try:
            return self.__polynomial
        except  AttributeError:
            R = polynomial_ring.PolynomialRing(FiniteField(self.characteristic()), name)
            self.__polynomial = R(self._pari_modulus())
        return self.__polynomial

    def __hash__(self):
        """
        Return the hash of this field.

        EXAMPLES:
            sage: hash(GF(3,'b'))
            904200654
            sage: hash(GF(3,'a'))
            904200654
            sage: hash(GF(9,'a'))
            -443918504
            sage: hash(GF(9,'b'))
            419125555
        """
        return hash((self.__order, self.variable_name(), self.__modulus))

class FiniteField_prime_modn(FiniteField_generic, integer_mod_ring.IntegerModRing_generic):
    def __init__(self, p, name=None):
        p = integer.Integer(p)
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        integer_mod_ring.IntegerModRing_generic.__init__(self, p)
        self.__char = p
        self.__gen = self(1)  # self(int(pari.pari(p).znprimroot().lift()))

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
            if K == self:
                return x
        raise TypeError, "no canonical coercion of x"

    def characteristic(self):
        return self.__char

    def modulus(self):
        try:
            return self.__modulus
        except AttributeError:
            x = polynomial_ring.PolynomialRing(self, 'x').gen()
            self.__modulus = x - 1
        return self.__modulus

    def is_prime_field(self):
        return True

    def is_prime(self):
        return True

    def polynomial(self, name):
        try:
            return self.__polynomial
        except  AttributeError:
            self.__polynomial = polynomial_ring.PolynomialRing(self, name)([0,1])
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
            sage: FiniteField(3**20, 'a').degree()
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
        RuntimeError: requested conway polynomial not in database.
    """
    (p,n)=(int(p),int(n))
    R = polynomial_ring.PolynomialRing(GF(p), 'x')
    try:
        return R(sage.databases.conway.ConwayPolynomials()[p][n])
    except KeyError:
        raise RuntimeError, "requested conway polynomial not in database."

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
        sage: F = GF(13, 'a')
        sage: F(x)
        2
        sage: F(gap('0*Z(13)'))
        0
        sage: F = GF(13^2, 'a')
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
        g = int(sage.interfaces.gap.gap.eval('Int(Z(%s))'%q))
    else:
        g = K.multiplicative_generator()
    return F(K(g**e))

