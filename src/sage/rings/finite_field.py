"""
Finite Fields


EXAMPLES:
Finite Fields support iteration, starting with 0.
    sage: k = GF(9, 'a')
    sage: for i,x in enumerate(k):  print i,x
    0 0
    1 2*a
    2 a + 1
    3 a + 2
    4 2
    5 a
    6 2*a + 2
    7 2*a + 1
    8 1
    sage: for a in GF(5):
    ...    print a
    0
    1
    2
    3
    4

We output the base rings of several finite fields.
    sage: k = GF(3); type(k)
    <class 'sage.rings.finite_field.FiniteField_prime_modn'>
    sage: k.base_ring()
    Finite Field of size 3

    sage: k = GF(9,'alpha'); type(k)
    <type 'sage.rings.finite_field_givaro.FiniteField_givaro'>
    sage: k.base_ring()
    Finite Field of size 3

    sage: k = GF(3^40,'b'); type(k)
    <class 'sage.rings.finite_field_ext_pari.FiniteField_ext_pari'>
    sage: k.base_ring()
    Finite Field of size 3

Further examples:
    sage: GF(2).is_field()
    True
    sage: GF(next_prime(10^20)).is_field()
    True
    sage: GF(19^20,'a').is_field()
    True
    sage: GF(8,'a').is_field()
    True
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
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

import random
import weakref

import arith

import integer
import rational
import integer_mod

import integer_mod_ring
from ring import is_FiniteField
from ring import FiniteField as FiniteField_generic
from finite_field_givaro import FiniteField_givaro

import polynomial.polynomial_ring as polynomial_ring
import polynomial.polynomial_element as polynomial_element
import polynomial.multi_polynomial_element as multi_polynomial_element

from sage.structure.parent_gens import normalize_names, ParentWithGens

import sage.interfaces.gap
import sage.databases.conway

cache = {}

def FiniteField(order, name=None, modulus=None, names=None,
                elem_cache=False, check_irreducible=True, *args, **kwds):
    """
    Return the globally unique finite field of given order with generator
    labeled by the given name and possibly with given modulus.

    INPUT:
        order --   int
        name --    string; must be specified if not a prime field
        modulus -- (optional) defining polynomial for field, i.e.,
                   generator of the field will be a root of this
                   polynomial; if not specified the choice of
                   definining polynomials can be arbitrary.
        elem_cache -- cache all elements to avoid creation time  (default: order<500)
        check_irreducible -- verify that the polynomial modulus is irreducible
        args -- additional parameters passed to finite field implementations
        kwds -- additional keyword parameters passed to finite field implementations

    ALIAS:
        You can also use GF instead of FiniteField -- they are identical.

    EXAMPLES:
        sage: k.<a> = FiniteField(9); k
        Finite Field in a of size 3^2
        sage: parent(a)
        Finite Field in a of size 3^2
        sage: charpoly(a, 'y')
        y^2 + 2*y + 2

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x +1 )
        sage: f = K.modulus(); f
        x^5 + 4*x + 1
        sage: type(f)
        <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_mod_p'>

    The modulus must be irreducible:
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x )
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not

    You can't accidently fool the constructor into thinking the
    modulus is irreducible when it isn't mod p, since it actually
    tests irreducibility modulo p.

        sage: F.<x> = QQ[]
        sage: factor(x^5+2)
        x^5 + 2
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 + 2 )
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not

    If you wish to live dangerously, you can tell the constructor not
    to test irreducibility using check_irreducible=False, but this
    can easily lead to crashes and hangs -- so do not do it unless
    you know that the modulus really is irreducible!

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**2, name='a', modulus=x^2 + 2, check_irreducible=False)

    For example, you may print finite field elements as integers via
    the Givaro implementation. But the constructor parameter to allow
    this is not passed to the actual implementation so far.

        sage: k.<a> = GF(2^8,repr='int')
        sage: a
        2

    The order of a finite field must be a prime power:
        sage: GF(100)
        Traceback (most recent call last):
        ...
        ValueError: order of finite field must be a prime power
    """
    if not names is None: name = names
    order = int(order)
    name = normalize_names(1,name)

    if elem_cache is None:
        elem_cache = order < 500

    key = (order, name, modulus, str([args, kwds]))
    if cache.has_key(key):
        K = cache[key]()
        if not K is None:
            return K
    if arith.is_prime(order):
        K = FiniteField_prime_modn(order,*args,**kwds)
    else:
        if not arith.is_prime_power(order):
            raise ValueError, "order of finite field must be a prime power"
        if check_irreducible and polynomial_element.is_Polynomial(modulus):
            if modulus.parent().base_ring().characteristic() == 0:
                p = arith.factor(order)[0][0]
                modulus = modulus.change_ring(FiniteField(p))
            if not modulus.is_irreducible():
                raise ValueError, "finite field modulus must be irreducible but it is not"
        if name is None:
            raise TypeError, "you must specify the generator name"
        if order < zech_log_bound:
            # DO *NOT* use for prime subfield, since that would lead to
            # a circular reference in the call to ParentWithGens in the
            # __init__ method.
            K = FiniteField_givaro(order, name, modulus, cache=elem_cache, *args,**kwds)
        else:
            if integer.Integer(order).factor()[0][0] == 2:
                from finite_field_ntl_gf2e import FiniteField_ntl_gf2e
                K = FiniteField_ntl_gf2e(order, name, modulus, *args, **kwds)
            else:
                from finite_field_ext_pari import FiniteField_ext_pari
                K = FiniteField_ext_pari(order, name, modulus, *args, **kwds)

    cache[key] = weakref.ref(K)
    return K


def is_PrimeFiniteField(x):
    """
    Returns True if x is a prime finite field (which is a specific
    data type).

    EXAMPLES:
        sage: is_PrimeFiniteField(QQ)
        False
        sage: is_PrimeFiniteField(GF(7))
        True
        sage: is_PrimeFiniteField(GF(7^2,'a'))
        False
        sage: is_PrimeFiniteField(GF(next_prime(10^90,proof=False)))
        True
    """
    return isinstance(x, FiniteField_prime_modn)

GF = FiniteField

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
        K = FiniteField(q, F.variable_name())
    if s.find(')^') == -1:
        e = 1
    else:
        e = int(s[i2+2:])
    if F.degree() == 1:
        g = int(sage.interfaces.gap.gap.eval('Int(Z(%s))'%q))
    else:
        g = K.multiplicative_generator()
    return F(K(g**e))


class FiniteField_prime_modn(FiniteField_generic, integer_mod_ring.IntegerModRing_generic):
    def __init__(self, p, name=None):
        p = integer.Integer(p)
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        integer_mod_ring.IntegerModRing_generic.__init__(self, p)
        self._kwargs = {}
        self.__char = p
        self.__gen = self(1)  # self(int(pari.pari(p).znprimroot().lift()))
        ParentWithGens.__init__(self, self, ('x',), normalize=False)

    def __cmp__(self, other):
        if not isinstance(other, FiniteField_prime_modn):
            return cmp(type(self), type(other))
        return cmp(self.__char, other.__char)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        This is called implicitly by the hom constructor.

        EXAMPLES:
            sage: k = GF(73^2,'a')
            sage: f = k.modulus()
            sage: r = f.change_ring(k).roots()
            sage: k.hom([r[0][0]])
            Ring endomorphism of Finite Field in a of size 73^2
              Defn: a |--> 72*a + 3
        """
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, integer.Integer)):
            return self(x)
        if isinstance(x, integer_mod.IntegerMod_abstract) and \
               x.parent().characteristic() == self.characteristic():
            return self(x)
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

    def polynomial(self, name=None):
        if name is None:
            name = self.variable_name()
        try:
            return self.__polynomial[name]
        except  AttributeError:
            R = polynomial_ring.PolynomialRing(FiniteField(self.characteristic()), name)
            f = polynomial_ring.PolynomialRing(self, name)([0,1])
            try:
                self.__polynomial[name] = f
            except (KeyError, AttributeError):
                self.__polynomial = {}
                self.__polynomial[name] = f
            return f

    def order(self):
        return self.__char

    def gen(self, n=0):
        """
        Return generator of this finite field.

        EXAMPLES:
            sage: k = GF(13)
            sage: k.gen()
            1
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: only one generator
        """
        if n != 0:
            raise IndexError, "only one generator"
        return self.__gen

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
            sage: FiniteField(3^20, 'a').degree()
            20
        """
        return 1

zech_log_bound = 2**16
