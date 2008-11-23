r"""
Number Fields

AUTHORS:
   -- William Stein (2004, 2005): initial version
   -- Steven Sivek (2006-05-12): added support for relative extensions
   -- William Stein (2007-09-04): major rewrite and documentation

NOTE:

    Unlike in PARI/GP, class group computations *in SAGE* do *not* by
    default assume the Generalized Riemann Hypothesis.  To do class
    groups computations not provably correctly you must often pass the
    flag proof=False to functions or call the function
    \code{proof.number_field(False)}.  It can easily take 1000's of
    times longer to do computations with \code{proof=True} (the
    default).

This example follows one in the Magma reference manual:
    sage: K.<y> = NumberField(x^4 - 420*x^2 + 40000)
    sage: z = y^5/11; z
    420/11*y^3 - 40000/11*y
    sage: R.<y> = PolynomialRing(K)
    sage: f = y^2 + y + 1
    sage: L.<a> = K.extension(f); L
    Number Field in a with defining polynomial y^2 + y + 1 over its base field
    sage: KL.<b> = NumberField([x^4 - 420*x^2 + 40000, x^2 + x + 1]); KL
    Number Field in b0 with defining polynomial x^4 - 420*x^2 + 40000 over its base field

We do some arithmetic in a tower of relative number fields:
    sage: K.<cuberoot2> = NumberField(x^3 - 2)
    sage: L.<cuberoot3> = K.extension(x^3 - 3)
    sage: S.<sqrt2> = L.extension(x^2 - 2)
    sage: S
    Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field
    sage: sqrt2 * cuberoot3
    cuberoot3*sqrt2
    sage: (sqrt2 + cuberoot3)^5
    (20*cuberoot3^2 + 15*cuberoot3 + 4)*sqrt2 + 3*cuberoot3^2 + 20*cuberoot3 + 60
    sage: cuberoot2 + cuberoot3
    cuberoot3 + cuberoot2
    sage: cuberoot2 + cuberoot3 + sqrt2
    sqrt2 + cuberoot3 + cuberoot2
    sage: (cuberoot2 + cuberoot3 + sqrt2)^2
    (2*cuberoot3 + 2*cuberoot2)*sqrt2 + cuberoot3^2 + 2*cuberoot2*cuberoot3 + cuberoot2^2 + 2
    sage: cuberoot2 + sqrt2
    sqrt2 + cuberoot2
    sage: a = S(cuberoot2); a
    cuberoot2
    sage: a.parent()
    Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field

WARNING: Doing arithmetic in towers of relative fields that depends on
canonical coercions is currently VERY SLOW.  It is much better to
explicitly coerce all elements into a common field, then do arithmetic
with them there (which is quite fast).

TESTS:
    sage: y = polygen(QQ,'y'); K.<beta> = NumberField([y^3 - 3, y^2 - 2])
    sage: K(y^10)
    (-3024*beta1 + 1530)*beta0^2 + (-2320*beta1 + 5067)*beta0 - 3150*beta1 + 7592
"""


# TODO:
#   * relative over relative

#*****************************************************************************
#       Copyright (C) 2004, 2005, 2006, 2007 William Stein <wstein@gmail.com>
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

from __future__ import with_statement
from sage.structure.parent_gens import localvars

# There will be one running instance of GP for all
# number field calculations that use the interpreter.
from sage.interfaces.gp import Gp

import sage.libs.ntl.all as ntl
import sage.libs.pari.all as pari
import sage.interfaces.gap
import sage.misc.preparser
import sage.rings.arith

import sage.rings.complex_field
import sage.rings.real_mpfr
import sage.rings.real_mpfi
import sage.rings.complex_double
import sage.rings.real_double

import sage.rings.ring
from sage.misc.latex import latex_variable_name, latex_varify

from class_group import ClassGroup
from galois_group import GaloisGroup
#import order

from sage.structure.element import is_Element
from sage.structure.sequence import Sequence

import sage.structure.parent_gens

from sage.structure.proof.proof import get_flag
import maps
from itertools import count, izip

from sage.rings.integer_ring import IntegerRing

def proof_flag(t):
    """
    Used for easily determining the correct proof flag to use.

    Returns t if t is not None, otherwise returns the system-wide
    proof-flag for number fields (default: True).

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import proof_flag
        sage: proof_flag(True)
        True
        sage: proof_flag(False)
        False
        sage: proof_flag(None)
        True
        sage: proof_flag("banana")
        'banana'
    """
    return get_flag(t, "number_field")

_gp = None
def gp():
    """
    Return the unique copy of the gp (PARI) interpreter
    used for number field computations.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import gp
        sage: gp()
        GP/PARI interpreter
    """
    global _gp
    if not _gp is None:
        return _gp
    else:
        _gp = Gp()
        return _gp

import operator

import weakref

from sage.misc.latex import latex

import sage.rings.arith as arith
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.infinity as infinity
import sage.rings.rational as rational
import sage.rings.integer as integer
import sage.rings.polynomial.polynomial_ring as polynomial_ring
import sage.rings.polynomial.polynomial_element as polynomial_element
import sage.rings.ideal as ideal
import sage.rings.complex_field
import sage.groups.abelian_gps.abelian_group
import sage.rings.complex_interval_field

from sage.structure.parent_gens import ParentWithGens
import number_field_element
import number_field_element_quadratic
from number_field_ideal import convert_from_zk_basis, NumberFieldIdeal, is_NumberFieldIdeal, NumberFieldFractionalIdeal
from sage.rings.number_field.number_field_ideal_rel import NumberFieldFractionalIdeal_rel
from sage.libs.all import pari, pari_gen

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()
RIF = sage.rings.real_mpfi.RealIntervalField()
CIF = sage.rings.complex_interval_field.ComplexIntervalField()

_nf_cache = {}
def NumberField(polynomial, name=None, check=True, names=None, cache=True):
    r"""
    Return \emph{the} number field defined by the given irreducible
    polynomial and with variable with the given name.  If check is
    True (the default), also verify that the defining polynomial is
    irreducible and over Q.

    INPUT:
        polynomial -- a polynomial over QQ or a number field, or
                      a list of polynomials.
        name -- a string (default: 'a'), the name of the generator
        check -- bool (default: True); do type checking and
                 irreducibility checking.

    EXAMPLES:
        sage: z = QQ['z'].0
        sage: K = NumberField(z^2 - 2,'s'); K
        Number Field in s with defining polynomial z^2 - 2
        sage: s = K.0; s
        s
        sage: s*s
        2
        sage: s^2
        2

    Constructing a relative number field
        sage: K.<a> = NumberField(x^2 - 2)
        sage: R.<t> = K[]
        sage: L.<b> = K.extension(t^3+t+a); L
        Number Field in b with defining polynomial t^3 + t + a over its base field
        sage: L.absolute_field('c')
        Number Field in c with defining polynomial x^6 + 2*x^4 + x^2 - 2
        sage: a*b
        a*b
        sage: L(a)
        a
        sage: L.lift_to_base(b^3 + b)
        -a

    Constructing another number field:
        sage: k.<i> = NumberField(x^2 + 1)
        sage: R.<z> = k[]
        sage: m.<j> = NumberField(z^3 + i*z + 3)
        sage: m
        Number Field in j with defining polynomial z^3 + i*z + 3 over its base field

    Number fields are globally unique.
        sage: K.<a>= NumberField(x^3-5)
        sage: a^3
        5
        sage: L.<a>= NumberField(x^3-5)
        sage: K is L
        True

    Having different defining polynomials makes them fields different:
        sage: x = polygen(QQ, 'x'); y = polygen(QQ, 'y')
        sage: k.<a> = NumberField(x^2 + 3)
        sage: m.<a> = NumberField(y^2 + 3)
        sage: k
        Number Field in a with defining polynomial x^2 + 3
        sage: m
        Number Field in a with defining polynomial y^2 + 3

    An example involving a variable name that defines a function in
    PARI:
        sage: theta = polygen(QQ, 'theta')
        sage: M.<z> = NumberField([theta^3 + 4, theta^2 + 3]); M
        Number Field in z0 with defining polynomial theta^3 + 4 over its base field

    TESTS:
        sage: x = QQ['x'].gen()
        sage: y = ZZ['y'].gen()
        sage: K = NumberField(x^3 + x + 3, 'a'); K
        Number Field in a with defining polynomial x^3 + x + 3
        sage: K.defining_polynomial().parent()
        Univariate Polynomial Ring in x over Rational Field

        sage: L = NumberField(y^3 + y + 3, 'a'); L
        Number Field in a with defining polynomial y^3 + y + 3
        sage: L.defining_polynomial().parent()
        Univariate Polynomial Ring in y over Rational Field
    """
    if name is None and names is None:
        raise TypeError, "You must specify the name of the generator."
    if not names is None:
        name = names

    if isinstance(polynomial, (list, tuple)):
        return NumberFieldTower(polynomial, name)

    name = sage.structure.parent_gens.normalize_names(1, name)

    if not isinstance(polynomial, polynomial_element.Polynomial):
        try:
            polynomial = polynomial.polynomial(QQ)
        except (AttributeError, TypeError):
            raise TypeError, "polynomial (=%s) must be a polynomial."%repr(polynomial)

    # convert ZZ to QQ
    R = polynomial.base_ring()
    Q = polynomial.parent().base_extend(R.fraction_field())
    polynomial = Q(polynomial)

    if cache:
        key = (polynomial, name)
        if _nf_cache.has_key(key):
            K = _nf_cache[key]()
            if not K is None: return K

    if isinstance(R, NumberField_generic):
        S = R.extension(polynomial, name, check=check)
        if cache:
            _nf_cache[key] = weakref.ref(S)
        return S

    if polynomial.degree() == 2:
        K = NumberField_quadratic(polynomial, name, check)
    else:
        K = NumberField_absolute(polynomial, name, None, check)

    if cache:
        _nf_cache[key] = weakref.ref(K)
    return K


def NumberFieldTower(v, names, check=True):
    """
    Return the tower of number fields defined by the polynomials or
    number fields in the list v.

    This is the field constructed first from v[0], then over that
    field from v[1], etc.  If all is False, then each v[i] must be
    irreducible over the previous fields.  Otherwise a list of all
    possible fields defined by all those polynomials is output.

    If names defines a variable name a, say, then the generators of
    the intermediate number fields are a0, a1, a2, ...

    INPUT:
        v -- a list of polynomials or number fields
        names -- variable names
        check -- bool (default: True) only relevant if all is False.
               Then check irreducibility of each input polynomial.

    OUTPUT:
        a single number field or a list of number fields

    EXAMPLES:
        sage: k.<a,b,c> = NumberField([x^2 + 1, x^2 + 3, x^2 + 5]); k
        Number Field in a with defining polynomial x^2 + 1 over its base field
        sage: a^2
        -1
        sage: b^2
        -3
        sage: c^2
        -5
        sage: (a+b+c)^2
        (2*b + 2*c)*a + 2*c*b - 9

    The Galois group is a product of 3 groups of order 2:
        sage: k.galois_group()
        Galois group PARI group [8, 1, 3, "E(8)=2[x]2[x]2"] of degree 8 of the Number Field in a with defining polynomial x^2 + 1 over its base field


    Repeatedly calling base_field allows us to descend the internally
    constructed tower of fields:
        sage: k.base_field()
        Number Field in b with defining polynomial x^2 + 3 over its base field
        sage: k.base_field().base_field()
        Number Field in c with defining polynomial x^2 + 5
        sage: k.base_field().base_field().base_field()
        Rational Field

    In the following example the second polynomial is reducible over the
    first, so we get an error:
        sage: v = NumberField([x^3 - 2, x^3 - 2], names='a')
        Traceback (most recent call last):
        ...
        ValueError: defining polynomial (x^3 - 2) must be irreducible

    We mix polynomial parent rings:
        sage: k.<y> = QQ[]
        sage: m = NumberField([y^3 - 3, x^2 + x + 1, y^3 + 2], 'beta')
        sage: m
        Number Field in beta0 with defining polynomial y^3 - 3 over its base field
        sage: m.base_field ()
        Number Field in beta1 with defining polynomial x^2 + x + 1 over its base field

    A tower of quadratic fields:
        sage: K.<a> = NumberField([x^2 + 3, x^2 + 2, x^2 + 1])
        sage: K
        Number Field in a0 with defining polynomial x^2 + 3 over its base field
        sage: K.base_field()
        Number Field in a1 with defining polynomial x^2 + 2 over its base field
        sage: K.base_field().base_field()
        Number Field in a2 with defining polynomial x^2 + 1

    A bigger tower of quadratic fields.
        sage: K.<a2,a3,a5,a7> = NumberField([x^2 + p for p in [2,3,5,7]]); K
        Number Field in a2 with defining polynomial x^2 + 2 over its base field
        sage: a2^2
        -2
        sage: a3^2
        -3
        sage: (a2+a3+a5+a7)^3
        ((6*a5 + 6*a7)*a3 + 6*a7*a5 - 47)*a2 + (6*a7*a5 - 45)*a3 - 41*a5 - 37*a7
    """
    # there is an "all" option below -- we do not use it, since
    # I couldn't get it to work with PARI reliably, and it isn't
    # very useful in practice.
    try:
        names = sage.structure.parent_gens.normalize_names(len(v), names)
    except IndexError:
        names = sage.structure.parent_gens.normalize_names(1, names)
        if len(v) > 1:
            names = ['%s%s'%(names[0], i) for i in range(len(v))]

    if not isinstance(v, (list, tuple)):
        raise TypeError, "v must be a list or tuple"
    if len(v) == 0:
        return QQ
    if len(v) == 1:
        #if all:
        #    f = v[0]
        #    if not isinstance(f, polynomial_element.Polynomial):
        #        f = QQ['x'](v[0])
        #    F = f.factor()
        #    return [NumberField(F[i][0], names='%s%s'%(names[0],str(i))) for i in range(len(F))]
        #else:
        return NumberField(v[0], names=names)
    f = v[0]
    w = NumberFieldTower(v[1:], names=names[1:])
    if isinstance(f, polynomial_element.Polynomial):
        var = f.name()
    else:
        var = 'x'

    name = names[0]
    #if all:
    #    R = w[0][var]  # polynomial ring
    #else:
    R = w[var]  # polynomial ring

    f = R(f)
    i = 0

    sep = chr(ord(name[0]) + 1)
    #if all:
    #    z = []
    #    for k in w:
    #        for g, _ in f.factor():
    #            z.append(k.extension(g, names='%s%s%s'%(name, sep, str(i)), check=False))
    #            i += 1
    #    return z
    #else:
    return w.extension(f, name, check=check)


def QuadraticField(D, names, check=True):
    """
    Return a quadratic field obtained by adjoining a square root of
    $D$ to the rational numbers, where $D$ is not a perfect square.

    INPUT:
        D -- a rational number
        name -- variable name
        check -- bool (default: True)

    OUTPUT:
        A number field defined by a quadratic polynomial.

    EXAMPLES:
        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3
        sage: K.<theta> = QuadraticField(3); K
        Number Field in theta with defining polynomial x^2 - 3
        sage: QuadraticField(9, 'a')
        Traceback (most recent call last):
        ...
        ValueError: D must not be a perfect square.
        sage: QuadraticField(9, 'a', check=False)
        Number Field in a with defining polynomial x^2 - 9

    Quadratic number fields derive from general number fields.
        sage: from sage.rings.number_field.number_field import is_NumberField
        sage: type(K)
        <class 'sage.rings.number_field.number_field.NumberField_quadratic'>
        sage: is_NumberField(K)
        True
    """
    D = QQ(D)
    if check:
        if D.is_square():
            raise ValueError, "D must not be a perfect square."
    R = QQ['x']
    f = R([-D, 0, 1])
    return NumberField(f, names, check=False)

def is_AbsoluteNumberField(x):
    """
    Return True if x is an absolute number field.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import is_AbsoluteNumberField
        sage: is_AbsoluteNumberField(NumberField(x^2+1,'a'))
        True
        sage: is_AbsoluteNumberField(NumberField([x^3 + 17, x^2+1],'a'))
        False

    The rationals are a number field, but they're not of the absolute number field class.
        sage: is_AbsoluteNumberField(QQ)
        False
    """
    return isinstance(x, NumberField_absolute)

def is_QuadraticField(x):
    r"""
    Return True if x is of the quadratic \emph{number} field type.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import is_QuadraticField
        sage: is_QuadraticField(QuadraticField(5,'a'))
        True
        sage: is_QuadraticField(NumberField(x^2 - 5, 'b'))
        True
        sage: is_QuadraticField(NumberField(x^3 - 5, 'b'))
        False

    A quadratic field specially refers to a number field, not a finite
    field:
        sage: is_QuadraticField(GF(9,'a'))
        False
    """
    return isinstance(x, NumberField_quadratic)

def is_RelativeNumberField(x):
    """
    Return True if x is a relative number field.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import is_RelativeNumberField
        sage: is_RelativeNumberField(NumberField(x^2+1,'a'))
        False
        sage: k.<a> = NumberField(x^3 - 2)
        sage: l.<b> = k.extension(x^3 - 3); l
        Number Field in b with defining polynomial x^3 - 3 over its base field
        sage: is_RelativeNumberField(l)
        True
        sage: is_RelativeNumberField(QQ)
        False
    """
    return isinstance(x, NumberField_relative)

_cyclo_cache = {}
def CyclotomicField(n, names=None):
    r"""
    Return the n-th cyclotomic field, where n is a positive integer.

    INPUT:
        n -- a positive integer
        names -- name of generator (optional -- defaults to zetan).

    EXAMPLES:
    We create the $7$th cyclotomic field $\QQ(\zeta_7)$ with the
    default generator name.
        sage: k = CyclotomicField(7); k
        Cyclotomic Field of order 7 and degree 6
        sage: k.gen()
        zeta7

    Cyclotomic fields are of a special type.
        sage: type(k)
        <class 'sage.rings.number_field.number_field.NumberField_cyclotomic'>

    We can specify a different generator name as follows.
        sage: k.<z7> = CyclotomicField(7); k
        Cyclotomic Field of order 7 and degree 6
        sage: k.gen()
        z7

    The $n$ must be an integer.
        sage: CyclotomicField(3/2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    The degree must be positive.
        sage: CyclotomicField(0)
        Traceback (most recent call last):
        ...
        ValueError: n (=0) must be a positive integer

    The special case $n=1$ does \emph{not} return the rational numbers:
        sage: CyclotomicField(1)
        Cyclotomic Field of order 1 and degree 1

        sage: cf6 = CyclotomicField(6) ; z6 = cf6.0
        sage: cf3 = CyclotomicField(3) ; z3 = cf3.0
        sage: cf3(z6)
        zeta3 + 1
        sage: cf6(z3)
        zeta6 - 1
        sage: cf9 = CyclotomicField(9) ; z9 = cf9.0
        sage: cf18 = CyclotomicField(18) ; z18 = cf18.0
        sage: cf18(z9)
        zeta18^2
        sage: cf9(z18)
        -zeta9^5
        sage: cf18(z3)
        zeta18^3 - 1
        sage: cf18(z6)
        zeta18^3
        sage: cf18(z6)**2
        zeta18^3 - 1
        sage: cf9(z3)
        zeta9^3
    """
    n = ZZ(n)
    if n <= 0:
        raise ValueError, "n (=%s) must be a positive integer"%n

    if names is None:
        names = "zeta%s"%n
    names = sage.structure.parent_gens.normalize_names(1, names)
    key = (n, names)
    if _cyclo_cache.has_key(key):
        K = _cyclo_cache[key]()
        if not K is None: return K
    K = NumberField_cyclotomic(n, names)
    _cyclo_cache[key] = weakref.ref(K)
    return K

def is_CyclotomicField(x):
    """
    Return True if x is a cyclotomic field, i.e., of the special
    cyclotomic field class.  This function does not return True
    for a number field that just happens to be isomorphic to a
    cyclotomic field.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import is_CyclotomicField
        sage: is_CyclotomicField(NumberField(x^2 + 1,'zeta4'))
        False
        sage: is_CyclotomicField(CyclotomicField(4))
        True
        sage: is_CyclotomicField(CyclotomicField(1))
        True
        sage: is_CyclotomicField(QQ)
        False
        sage: is_CyclotomicField(7)
        False
    """
    return isinstance(x, NumberField_cyclotomic)

import number_field_base

is_NumberField = number_field_base.is_NumberField

class NumberField_generic(number_field_base.NumberField):
    """
    EXAMPLES:
        sage: K.<a> = NumberField(x^3 - 2); K
        Number Field in a with defining polynomial x^3 - 2
        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, polynomial, name,
                 latex_name=None, check=True):
        """
        Create a number field.

        EXAMPLES:
            sage: NumberField(x^97 - 19, 'a')
            Number Field in a with defining polynomial x^97 - 19

        The defining polynomial must be irreducible:
            sage: K.<a> = NumberField(x^2 - 1)
            Traceback (most recent call last):
            ...
            ValueError: defining polynomial (x^2 - 1) must be irreducible

        If you use check=False, you avoid checking irreducibility of
        the defining polynomial, which can save time.
            sage: K.<a> = NumberField(x^2 - 1, check=False)

        It can also be dangerous:
            sage: (a-1)*(a+1)
            0

        TESTS:
            sage: NumberField(ZZ['x'].0^4 + 23, 'a')
            Number Field in a with defining polynomial x^4 + 23
            sage: NumberField(QQ['x'].0^4 + 23, 'a')
            Number Field in a with defining polynomial x^4 + 23
            sage: NumberField(GF(7)['x'].0^4 + 23, 'a')
            Traceback (most recent call last):
            ...
            TypeError: polynomial must be defined over rational field
        """
        ParentWithGens.__init__(self, QQ, name)
        if not isinstance(polynomial, polynomial_element.Polynomial):
            raise TypeError, "polynomial (=%s) must be a polynomial"%repr(polynomial)

        if check:
            if not polynomial.parent().base_ring() == QQ:
                raise TypeError, "polynomial must be defined over rational field"
            if not polynomial.is_monic():
                raise NotImplementedError, "number fields for non-monic polynomials not yet implemented."
            if not polynomial.is_irreducible():
                raise ValueError, "defining polynomial (%s) must be irreducible"%polynomial

        self._assign_names(name)
        if latex_name is None:
            self.__latex_variable_name = latex_variable_name(self.variable_name())
        else:
            self.__latex_variable_name = latex_name
        self.__polynomial = polynomial
        self.__pari_bnf_certified = False
        self._integral_basis_dict = {}

    def _Hom_(self, codomain, cat=None):
        """
        Return homset of homomorphisms from self to the number field codomain.

        The cat option is currently ignored.

        EXAMPLES:
        This function is implicitly called by the Hom method or function.
            sage: K.<i> = NumberField(x^2 + 1); K
            Number Field in i with defining polynomial x^2 + 1
            sage: K.Hom(K)
            Automorphism group of Number Field in i with defining polynomial x^2 + 1
            sage: Hom(K, QuadraticField(-1, 'b'))
            Set of field embeddings from Number Field in i with defining polynomial x^2 + 1 to Number Field in b with defining polynomial x^2 + 1
        """
        import morphism
        return morphism.NumberFieldHomset(self, codomain)

    def _set_structure(self, from_self, to_self, unsafe_force_change=False):
        """
        Internal function to set the structure fields of this number field.
        """
        # Note -- never call this on a cached number field, since
        # that could eventually lead to problems.
        if unsafe_force_change:
            self.__from_self = from_self
            self.__to_self = to_self
            return
        try:
            self.__from_self
        except AttributeError:
            self.__from_self = from_self
            self.__to_self = to_self
        else:
            raise ValueError, "number field structure is immutable."

    def structure(self):
        """
        Return fixed isomorphism or embedding structure on self.

        This is used to record various isomorphisms or embeddings
        that arise naturally in other constructions.

        EXAMPLES:
            sage: K.<z> = NumberField(x^2 + 3)
            sage: L.<a> = K.absolute_field(); L
            Number Field in a with defining polynomial x^2 + 3
            sage: L.structure()
            (Number field isomorphism from Number Field in a with defining polynomial x^2 + 3 to Number Field in z with defining polynomial x^2 + 3 given by variable name change,
             Number field isomorphism from Number Field in z with defining polynomial x^2 + 3 to Number Field in a with defining polynomial x^2 + 3 given by variable name change)
        """
        try:
            if self.__to_self is not None:
                return self.__from_self, self.__to_self
            else:
                return self.__from_self, self.__to_self
        except AttributeError:
            f = self.hom(self)
            self._set_structure(f,f)
            return f, f

    def primitive_element(self):
        r"""
        Return a primitive element for this field, i.e., an element
        that generates it over $\QQ$.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 2)
            sage: K.primitive_element()
            a
            sage: K.<a,b,c> = NumberField([x^2-2,x^2-3,x^2-5])
            sage: K.primitive_element()
            a - b + c
            sage: alpha = K.primitive_element(); alpha
            a - b + c
            sage: alpha.minpoly()
            x^2 + (2*b - 2*c)*x - 2*c*b + 6
            sage: alpha.absolute_minpoly()
            x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
        """
        try:
            return self.__primitive_element
        except AttributeError:
            pass
        K = self.absolute_field('a')
        from_K, to_K = K.structure()
        self.__primitive_element = from_K(K.gen())
        return self.__primitive_element

    def subfield(self, alpha, name=None):
        r"""
        Return an absolute number field K isomorphic to QQ(alpha) and
        a map from K to self that sends the generator of K to alpha.

        INPUT:
            alpha -- an element of self, or something that coerces
                     to an element of self.

        OUTPUT:
            K -- a number field
            from_K -- a homomorphism from K to self that sends the
                      generator of K to alpha.

        EXAMPLES:
            sage: K.<a> = NumberField(x^4 - 3); K
            Number Field in a with defining polynomial x^4 - 3
            sage: H, from_H = K.subfield(a^2, name='b')
            sage: H
            Number Field in b with defining polynomial x^2 - 3
            sage: from_H(H.0)
            a^2
            sage: from_H
            Ring morphism:
              From: Number Field in b with defining polynomial x^2 - 3
              To:   Number Field in a with defining polynomial x^4 - 3
              Defn: b |--> a^2

            sage: K.<z> = CyclotomicField(5)
            sage: K.subfield(z-z^2-z^3+z^4)
            (Number Field in z0 with defining polynomial x^2 - 5,
            Ring morphism:
            From: Number Field in z0 with defining polynomial x^2 - 5
            To:   Cyclotomic Field of order 5 and degree 4
            Defn: z0 |--> -2*z^3 - 2*z^2 - 1)

        You can also view a number field as having a different
        generator by just choosing the input to generate the
        whole field; for that it is better to use
        \code{self.change_generator}, which gives isomorphisms
        in both directions.
        """
        if name is None:
            name = self.variable_name() + '0'
        beta = self(alpha)
        f = beta.minpoly()
        K = NumberField(f, names=name)
        from_K = K.hom([beta])
        return K, from_K

    def change_generator(self, alpha, name=None):
        r"""
        Given the number field self, construct another isomorphic
        number field $K$ generated by the element alpha of self, along
        with isomorphisms from $K$ to self and from self to $K$.

        EXAMPLES:
            sage: K.<i> = NumberField(x^2 + 1); K
            Number Field in i with defining polynomial x^2 + 1
            sage: L.<i> = NumberField(x^2 + 1); L
            Number Field in i with defining polynomial x^2 + 1
            sage: K, from_K, to_K = L.change_generator(i/2 + 3)
            sage: K
            Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
            sage: from_K
            Ring morphism:
              From: Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
              To:   Number Field in i with defining polynomial x^2 + 1
              Defn: i0 |--> 1/2*i + 3
            sage: to_K
            Ring morphism:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Number Field in i0 with defining polynomial x^2 - 6*x + 37/4
              Defn: i |--> 2*i0 - 6

        We compute the image of the generator $\sqrt{-1}$ of $L$.
            sage: to_K(L.0)
            2*i0 - 6

        Note that he image is indeed a square root of -1.
            sage: to_K(L.0)^2
            -1
            sage: from_K(to_K(L.0))
            i
            sage: to_K(from_K(K.0))
            i0
        """
        alpha = self(alpha)
        K, from_K = self.subfield(alpha, name=name)
        if K.degree() != self.degree():
            raise ValueError, "alpha must generate a field of degree %s, but alpha generates a subfield of degree %s"%(self.degree(), K.degree())
        # Now compute to_K, which is an isomorphism
        # from self to K such that from_K(to_K(x)) == x for all x,
        # and to_K(from_K(y)) == y.
        # To do this, we must compute the image of self.gen()
        # under to_K.   This means writing self.gen() as a
        # polynomial in alpha, which is possible by the degree
        # check above.  This latter we do by linear algebra.
        phi = alpha.coordinates_in_terms_of_powers()
        c = phi(self.gen())
        to_K = self.hom([K(c)])
        return K, from_K, to_K

    def is_absolute(self):
        """
        Returns True if self is an absolute field.

        This function will be implemented in the derived classes.

        EXAMPLES:
            sage: K = CyclotomicField(5)
            sage: K.is_absolute()
            True
        """
        raise NotImplementedError

    def is_relative(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^10 - 2)
            sage: K.is_absolute()
            True
            sage: K.is_relative()
            False
        """
        return not self.is_absolute()

    def absolute_field(self, names):
        """
        Returns self as an absolute extension over QQ.

        OUTPUT:
            K -- this number field (since it is already absolute)

        Also, \code{K.structure()} returns from_K and to_K, where
        from_K is an isomorphism from K to self and to_K is an isomorphism
        from self to K.

        EXAMPLES:
            sage: K = CyclotomicField(5)
            sage: K.absolute_field('a')
            Number Field in a with defining polynomial x^4 + x^3 + x^2 + x + 1
        """
        try:
            return self.__absolute_field[names]
        except KeyError:
            pass
        except AttributeError:
            self.__absolute_field = {}
        K = NumberField(self.defining_polynomial(), names, cache=False)
        K._set_structure(maps.NameChangeMap(K, self), maps.NameChangeMap(self, K))
        self.__absolute_field[names] = K
        return K

    def is_isomorphic(self, other):
        """
        Return True if self is isomorphic as a number field to other.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 1)
            sage: m.<b> = NumberField(x^2 + 4)
            sage: k.is_isomorphic(m)
            True
            sage: m.<b> = NumberField(x^2 + 5)
            sage: k.is_isomorphic (m)
            False

            sage: k = NumberField(x^3 + 2, 'a')
            sage: k.is_isomorphic(NumberField((x+1/3)^3 + 2, 'b'))
            True
            sage: k.is_isomorphic(NumberField(x^3 + 4, 'b'))
            True
            sage: k.is_isomorphic(NumberField(x^3 + 5, 'b'))
            False
        """
        if not isinstance(other, NumberField_generic):
            raise ValueError, "other must be a generic number field."
        t = self.pari_polynomial().nfisisom(other.pari_polynomial())
        return t != 0

    def is_totally_real(self):
        """
        Return True if self is totally real, and False otherwise.

        Totally real means that every isomorphic embedding of self into the
        complex numbers has image contained in the real numbers.

        EXAMPLES:
            sage: NumberField(x^2+2, 'alpha').is_totally_real()
            False
            sage: NumberField(x^2-2, 'alpha').is_totally_real()
            True
            sage: NumberField(x^4-2, 'alpha').is_totally_real()
            False
        """
        return self.signature()[1] == 0

    def is_totally_imaginary(self):
        """
        Return True if self is totally imaginary, and False otherwise.

        Totally imaginary means that no isomorphic embedding of self into the
        complex numbers has image contained in the real numbers.

        EXAMPLES:
            sage: NumberField(x^2+2, 'alpha').is_totally_imaginary()
            True
            sage: NumberField(x^2-2, 'alpha').is_totally_imaginary()
            False
            sage: NumberField(x^4-2, 'alpha').is_totally_imaginary()
            False
        """
        return self.signature()[0] == 0

    def complex_embeddings(self, prec=53):
        r"""
        Return all homomorphisms of this number field into the
        approximate complex field with precision prec.

        If prec is 53 (the default), then the complex double field is
        used; otherwise the arbitrary precision (but slow) complex
        field is used.  If you want 53-bit arbitrary precision then
        do \code{self.embeddings(ComplexField(53))}.

        EXAMPLES:
            sage: k.<a> = NumberField(x^5 + x + 17)
            sage: v = k.complex_embeddings()
            sage: ls = [phi(k.0^2) for phi in v] ; ls # random order
            [2.97572074038...,
             -2.40889943716 + 1.90254105304*I,
             -2.40889943716 - 1.90254105304*I,
             0.921039066973 + 3.07553311885*I,
             0.921039066973 - 3.07553311885*I]
            sage: K.<a> = NumberField(x^3 + 2)
            sage: ls = K.complex_embeddings() ; ls # random order
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Complex Double Field
              Defn: a |--> -1.25992104989...,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Complex Double Field
              Defn: a |--> 0.629960524947 - 1.09112363597*I,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Complex Double Field
              Defn: a |--> 0.629960524947 + 1.09112363597*I
            ]
        """
        if prec == 53:
            CC = sage.rings.complex_double.CDF
        else:
            CC = sage.rings.complex_field.ComplexField(prec)
        return self.embeddings(CC)

    def real_embeddings(self, prec=53):
        r"""
        Return all homomorphisms of this number field into the
        approximate real field with precision prec.

        If prec is 53 (the default), then the real double field is
        used; otherwise the arbitrary precision (but slow) real field
        is used.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 2)
            sage: K.real_embeddings()
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Real Double Field
              Defn: a |--> -1.25992104989
            ]
             sage: K.real_embeddings(16)
             [
             Ring morphism:
               From: Number Field in a with defining polynomial x^3 + 2
               To:   Real Field with 16 bits of precision
               Defn: a |--> -1.260
             ]
            sage: K.real_embeddings(100)
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 + 2
              To:   Real Field with 100 bits of precision
              Defn: a |--> -1.2599210498948731647672106073
            ]
        """
        if prec == 53:
            K = sage.rings.real_double.RDF
        else:
            K = sage.rings.real_mpfr.RealField(prec)
        return self.embeddings(K)

    def latex_variable_name(self, name=None):
        """
        Return the latex representation of the variable name for
        this number field.

        EXAMPLES:
            sage: NumberField(x^2 + 3, 'a').latex_variable_name()
            'a'
            sage: NumberField(x^3 + 3, 'theta3').latex_variable_name()
            '\\theta_{3}'
            sage: CyclotomicField(5).latex_variable_name()
            '\\zeta_{5}'
        """
        if name is None:
            return self.__latex_variable_name
        else:
            self.__latex_variable_name = name

    def _repr_(self):
        """
        Return string representation of this number field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3)
            sage: k._repr_()
            'Number Field in a with defining polynomial x^13 - 2/3*x + 3'
        """
        return "Number Field in %s with defining polynomial %s"%(
                   self.variable_name(), self.polynomial())

    def _latex_(self):
        r"""
        Return latex representation of this number field.  This is viewed
        as a polynomial quotient ring over a field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3)
            sage: k._latex_()
            '\\mathbf{Q}[a]/(a^{13} - \\frac{2}{3} a + 3)'
            sage: latex(k)
            \mathbf{Q}[a]/(a^{13} - \frac{2}{3} a + 3)

        Numbered variables are often correctly typeset:
            sage: k.<theta25> = NumberField(x^25+x+1)
            sage: print k._latex_()
            \mathbf{Q}[\theta_{25}]/(\theta_{25}^{25} + \theta_{25} + 1)
        """
        return "%s[%s]/(%s)"%(latex(QQ), self.latex_variable_name(),
                              self.polynomial()._latex_(self.latex_variable_name()))

    def __call__(self, x):
        """
        Coerce x into this number field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 17)
            sage: K(a) is a
            True
            sage: K('a^2 + 2/3*a + 5')
            a^2 + 2/3*a + 5
            sage: K('1').parent()
            Number Field in a with defining polynomial x^3 + 17
            sage: K(3/5).parent()
            Number Field in a with defining polynomial x^3 + 17
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            K = x.parent()
            if K is self:
                return x
            elif K == self:
                return self._element_class(self, x.polynomial())
            elif isinstance(x, (number_field_element.OrderElement_absolute,
                                number_field_element.OrderElement_relative,
                                number_field_element_quadratic.OrderElement_quadratic)):
                L = K.number_field()
                if L is self:
                    return self._element_class(self, x)
                x = L(x)
            return self._coerce_from_other_number_field(x)
        elif isinstance(x,str):
            return self._coerce_from_str(x)
        elif isinstance(x, (sage.modules.vector_integer_dense.Vector_integer_dense,
                            sage.modules.vector_rational_dense.Vector_rational_dense)):
            if len(x) != self.degree():
                raise ValueError, "vector must be of length equal to the degree of this number field"
            return sum([ x[i]*self.gen(0)**i for i in range(self.degree()) ])
        return self._coerce_non_number_field_element_in(x)

    def _coerce_from_str(self, x):
        """
        Coerce a string representation of an element of this
        number field into this number field.

        INPUT:
            x -- string

        EXAMPLES:
            sage: k.<theta25> = NumberField(x^3+(2/3)*x+1)
            sage: k._coerce_from_str('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1

        This function is called by the coerce method when it gets a string
        as input:
            sage: k('theta25^3 + (1/3)*theta25')
            -1/3*theta25 - 1
        """
        # provide string coercion, as
        # for finite fields
        w = sage.misc.all.sage_eval(x,locals=\
                                  {self.variable_name():self.gen()})
        if not (is_Element(w) and w.parent() is self):
            return self(w)
        else:
            return w

    def _coerce_from_other_number_field(self, x):
        """
        Coerce a number field element x into this number field.

        In most cases this currently doesn't work (since it is
        barely implemented) -- it only works for constants.

        INPUT:
            x -- an element of some number field

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 2)
            sage: L.<b> = NumberField(x^2 + 1)
            sage: K._coerce_from_other_number_field(L(2/3))
            2/3
        """
        f = x.polynomial()
        if f.degree() <= 0:
            return self._element_class(self, f[0])
        # todo: more general coercion if embedding have been asserted
        raise TypeError, "Cannot coerce element into this number field"

    def _coerce_non_number_field_element_in(self, x):
        """
        Coerce a non-number field element x into this number field.

        INPUT:
            x -- a non number field element x, e.g., a list, integer,
            rational, or polynomial.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 2/3)
            sage: K._coerce_non_number_field_element_in(-7/8)
            -7/8
            sage: K._coerce_non_number_field_element_in([1,2,3])
            3*a^2 + 2*a + 1

        The list is just turned into a polynomial in the generator.
            sage: K._coerce_non_number_field_element_in([0,0,0,1,1])
            -2/3*a - 2/3

        Any polynomial whose coefficients can be coerced to rationals will
        coerce, e.g., this one in characteristic 7.
            sage: f = GF(7)['y']([1,2,3]); f
            3*y^2 + 2*y + 1
            sage: K._coerce_non_number_field_element_in(f)
            3*a^2 + 2*a + 1

       But not this one over a field of order 27.
            sage: F27.<g> = GF(27)
            sage: f = F27['z']([g^2, 2*g, 1]); f
            z^2 + 2*g*z + g^2
            sage: K._coerce_non_number_field_element_in(f)
            Traceback (most recent call last):
            ...
            TypeError
        """
        if isinstance(x, (int, long, rational.Rational,
                              integer.Integer, pari_gen,
                              list)):
            return self._element_class(self, x)

        try:
            if isinstance(x, polynomial_element.Polynomial):
                return self._element_class(self, x)

            return self._element_class(self, x._rational_())
        except (TypeError, AttributeError), msg:
            pass
        raise TypeError

    def _coerce_impl(self, x):
        """
        Canonical coercion of x into self.

        Currently integers, rationals, and this field itself coerce
        canonical into this field.

        EXAMPLES:
            sage: S.<y> = NumberField(x^3 + x + 1)
            sage: S._coerce_impl(int(4))
            4
            sage: S._coerce_impl(long(7))
            7
            sage: S._coerce_impl(-Integer(2))
            -2
            sage: z = S._coerce_impl(-7/8); z, type(z)
            (-7/8, <type 'sage.rings.number_field.number_field_element.NumberFieldElement_absolute'>)
            sage: S._coerce_impl(y) is y
            True

        There are situations for which one might imagine canonical
        coercion could make sense (at least after fixing choices), but
        which aren't yet implemented:
            sage: K.<a> = QuadraticField(2)
            sage: K._coerce_impl(sqrt(2))
            Traceback (most recent call last):
            ...
            TypeError
        """
        if isinstance(x, (rational.Rational, integer.Integer, int, long)):
            return self._element_class(self, x)
        elif isinstance(x, number_field_element.NumberFieldElement):
            from sage.rings.number_field.order import is_NumberFieldOrder
            if x.parent() is self:
                return x
            elif is_NumberFieldOrder(x.parent()) and x.parent().number_field() is self:
                return self._element_class(self, x.polynomial())
            elif x.parent() == self:
                return self._element_class(self, x.list())
        raise TypeError

    def category(self):
        """
        Return the category of number fields.

        EXAMPLES:
            sage: NumberField(x^2 + 3, 'a').category()
            Category of number fields
            sage: category(NumberField(x^2 + 3, 'a'))
            Category of number fields

        The special types of number fields, e.g., quadratic fields,
        don't have their own category:
            sage: QuadraticField(2,'d').category()
            Category of number fields
        """
        from sage.categories.all import NumberFields
        return NumberFields()

    def __cmp__(self, other):
        """
        Compare a number field with something else.

        INPUT:
            other -- arbitrary Python object.

        If other is not a number field, then the types of self and
        other are compared.  If both are number fields, then the
        variable names are compared.  If those are the same, then the
        underlying defining polynomials are compared.  If the
        polynomials are the same, the number fields are considered
        ``equal'', but need not be identical.  Coercion between equal
        number fields is allowed.

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 + 2); m.<b> = NumberField(x^3 + 2)
            sage: cmp(k,m)
            -1
            sage: cmp(m,k)
            1
            sage: k == QQ
            False
            sage: k.<a> = NumberField(x^3 + 2); m.<a> = NumberField(x^3 + 2)
            sage: k is m
            True
            sage: m = loads(dumps(k))
            sage: k is m
            False
            sage: k == m
            True

        TESTS:
            sage: x = QQ['x'].gen()
            sage: y = ZZ['y'].gen()
            sage: K = NumberField(x^3 + x + 3, 'a'); K
            Number Field in a with defining polynomial x^3 + x + 3
            sage: K.defining_polynomial().parent()
            Univariate Polynomial Ring in x over Rational Field

            sage: L = NumberField(y^3 + y + 3, 'a'); L
            Number Field in a with defining polynomial y^3 + y + 3
            sage: L.defining_polynomial().parent()
            Univariate Polynomial Ring in y over Rational Field
            sage: L == K
            True

            sage: NumberField(ZZ['x'].0^4 + 23, 'a') == NumberField(ZZ['y'].0^4 + 23, 'a')
            True
            sage: NumberField(ZZ['x'].0^4 + 23, 'a') == NumberField(QQ['y'].0^4 + 23, 'a')
            True
            sage: NumberField(QQ['x'].0^4 + 23, 'a') == NumberField(QQ['y'].0^4 + 23, 'a')
            True

            sage: x = var('x'); y = ZZ['y'].gen()
            sage: NumberField(x^3 + x + 5, 'a') == NumberField(y^3 + y + 5, 'a')
            True
            sage: NumberField(x^3 + x + 5, 'a') == NumberField(y^4 + y + 5, 'a')
            False
            sage: NumberField(x^3 + x + 5, 'a') == NumberField(x^3 + x + 5, 'b')
            False
        """
        if not isinstance(other, NumberField_generic):
            return cmp(type(self), type(other))
        c = cmp(self.variable_name(), other.variable_name())
        if c: return c
        # compare coefficients so that the polynomial variable does not count
        return cmp(list(self.__polynomial), list(other.__polynomial))

    def _ideal_class_(self):
        """
        Return the Python class used in defining the zero ideal of the
        ring of integers of this number field.

        This function is required by the general ring/ideal machinery.
        The value defined here is the default value for all number
        fields.

        EXAMPLES:
            sage: NumberField(x^2 + 2, 'c')._ideal_class_()
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldIdeal'>
        """
        return sage.rings.number_field.number_field_ideal.NumberFieldIdeal

    def _fractional_ideal_class_(self):
        """
        Return the Python class used in defining fractional ideals of
        the ring of integers of this number field.

        This function is required by the general ring/ideal machinery.
        The value defined here is the default value for all number
        fields *except* relative number fields; this function is
        overridden by one of the same name on class
        NumberField_relative.

        EXAMPLES:
            sage: NumberField(x^2 + 2, 'c')._fractional_ideal_class_()
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal'>
        """
        return sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal

    def ideal(self, *gens, **kwds):
        """
        K.ideal() returns a fractional ideal of the field, except for the
        zero ideal which is not a fractional ideal.

        EXAMPLES:
            sage: K.<i>=NumberField(x^2+1)
            sage: K.ideal(2)
            Fractional ideal (2)
            sage: K.ideal(2+i)
            Fractional ideal (i + 2)
            sage: K.ideal(0)
            Ideal (0) of Number Field in i with defining polynomial x^2 + 1
        """
        try:
            return self.fractional_ideal(*gens, **kwds)
        except ValueError:
            return sage.rings.ring.Ring.ideal(self, gens, **kwds)

    def fractional_ideal(self, *gens, **kwds):
        r"""
        Return the ideal in $\mathcal{O}_K$ generated by gens.  This
        overrides the \code{sage.rings.ring.Field} method to use the
        \code{sage.rings.ring.Ring} one instead, since we're not really
        concerned with ideals in a field but in its ring of integers.

        INPUT:
            gens -- a list of generators, or a number field ideal.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-2)
            sage: K.fractional_ideal([1/a])
            Fractional ideal (1/2*a^2)

        One can also input in a number field ideal itself.
            sage: K.fractional_ideal(K.ideal(a))
            Fractional ideal (a)

        The zero ideal is not a fractional ideal!
            sage: K.fractional_ideal(0)
            Traceback (most recent call last):
            ...
            ValueError: gens must have a nonzero element (zero ideal is not a fractional ideal)

        """
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        if len(gens) == 1 and isinstance(gens[0],self._fractional_ideal_class_()):
            I = gens[0]
            if I.number_field() is self:
                return I
            else:
                gens = I.gens()
        return self._fractional_ideal_class_()(self, gens, **kwds)

    def ideals_of_bdd_norm(self, bound):
        """
        All integral ideals of bounded norm.

        INPUT:
            bound -- a positive integer

        OUTPUT:
            A dict of all integral ideals I such that Norm(I) <= bound,
            keyed by norm.

        EXAMPLE:
            sage: K.<a> = NumberField(x^2 + 23)
            sage: d = K.ideals_of_bdd_norm(10)
            sage: for n in d:
            ...       print n
            ...       for I in d[n]:
            ...           print I
            1
            Fractional ideal (1)
            2
            Fractional ideal (2, 1/2*a - 1/2)
            Fractional ideal (2, 1/2*a + 1/2)
            3
            Fractional ideal (3, -1/2*a + 1/2)
            Fractional ideal (3, -1/2*a - 1/2)
            4
            Fractional ideal (4, 1/2*a + 3/2)
            Fractional ideal (2)
            Fractional ideal (4, 1/2*a + 5/2)
            5
            6
            Fractional ideal (-1/2*a + 1/2)
            Fractional ideal (6, 1/2*a + 5/2)
            Fractional ideal (6, 1/2*a + 7/2)
            Fractional ideal (1/2*a + 1/2)
            7
            8
            Fractional ideal (-1/2*a - 3/2)
            Fractional ideal (4, a - 1)
            Fractional ideal (4, a + 1)
            Fractional ideal (1/2*a - 3/2)
            9
            Fractional ideal (9, 1/2*a + 11/2)
            Fractional ideal (3)
            Fractional ideal (9, 1/2*a + 7/2)
            10

        """
        from sage.rings.number_field.number_field_ideal import convert_from_zk_basis
        hnf_ideals = pari('ideallist(%s, %d)'%(self.pari_nf(),bound))
        d = {}
        for i in xrange(bound):
            d[i+1] = [self.ideal([ self(generator) for generator in convert_from_zk_basis(self, hnf_I) ]) for hnf_I in hnf_ideals[i]]
        return d

    def primes_above(self, x, degree=None):
        r"""
        Return prime ideals of self lying over x.

        INPUT:
            -- x: usually an element or ideal of self.  It should be such that
               self.ideal(x) is sensible.  This excludes x=0.
            -- degree (default: None): None or an integer.  If None, find all
               primes above x of any degree.  If an integer, find all primes
               above x such that the resulting residue field has exactly this
               degree.

        OUTPUT:
            A list of prime ideals of self lying over x.  If degree is
            specified and no such ideal exists, returns the empty list.

        WARNING: at this time we factor the ideal x, which may not be
        supported for relative number fields.

        EXAMPLES:
            sage: x = ZZ['x'].gen()
            sage: F.<t> = NumberField(x^3 - 2)

            sage: P2s = F.primes_above(2)
            sage: P2s # random
            [Fractional ideal (-t)]
            sage: all(2 in P2 for P2 in P2s)
            True
            sage: all(P2.is_prime() for P2 in P2s)
            True
            sage: [ P2.norm() for P2 in P2s ]
            [2]

            sage: P3s = F.primes_above(3)
            sage: P3s # random
            [Fractional ideal (t + 1)]
            sage: all(3 in P3 for P3 in P3s)
            True
            sage: all(P3.is_prime() for P3 in P3s)
            True
            sage: [ P3.norm() for P3 in P3s ]
            [3]

            The ideal (3) is totally ramified in F, so there is no degree 2
            prime above 3:

            sage: F.primes_above(3, degree=2)
            []
            sage: [ id.residue_class_degree() for id, _ in F.ideal(3).factor() ]
            [1]

            Asking for a specific degree works:

            sage: P5_1s = F.primes_above(5, degree=1)
            sage: P5_1s # random
            [Fractional ideal (-t^2 - 1)]
            sage: P5_1 = P5_1s[0]; P5_1.residue_class_degree()
            1

            sage: P5_2s = F.primes_above(5, degree=2)
            sage: P5_2s # random
            [Fractional ideal (t^2 - 2*t - 1)]
            sage: P5_2 = P5_2s[0]; P5_2.residue_class_degree()
            2

        TESTS:
            It doesn't make sense to factor the ideal (0):

            sage: F.primes_above(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'factor'

            Sage can't factor ideals over extension fields yet:

            sage: G = F.extension(x^2 - 11, 'b')
            sage: G.primes_above(13)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if degree is not None:
            degree = ZZ(degree)
        ideal = self.ideal(x)
        facs = [ (id.residue_class_degree(), id) for id, _ in ideal.factor() ]
        facs.sort() # sorts on residue_class_degree(), lowest first
        if degree is None:
            return [ id for d, id in facs ]
        else:
            return [ id for d, id in facs if d == degree ]

    def prime_above(self, x, degree=None):
        r"""
        Return a prime ideal of self lying over x.

        INPUT:
            -- x: usually an element or ideal of self.  It should be such that
               self.ideal(x) is sensible.  This excludes x=0.
            -- degree (default: None): None or an integer.  If one, find a
               prime above x of any degree.  If an integer, find a prime above
               x such that the resulting residue field has exactly this degree.

        OUTPUT:
            A prime ideal of self lying over x.  If degree is specified and no
            such ideal exists, raises a ValueError.

        WARNING: at this time we factor the ideal x, which may not be
        supported for relative number fields.

        EXAMPLES:
            sage: x = ZZ['x'].gen()
            sage: F.<t> = NumberField(x^3 - 2)

            sage: P2 = F.prime_above(2)
            sage: P2 # random
            Fractional ideal (-t)
            sage: 2 in P2
            True
            sage: P2.is_prime()
            True
            sage: P2.norm()
            2

            sage: P3 = F.prime_above(3)
            sage: P3 # random
            Fractional ideal (t + 1)
            sage: 3 in P3
            True
            sage: P3.is_prime()
            True
            sage: P3.norm()
            3

            The ideal (3) is totally ramified in F, so there is no degree 2
            prime above 3:

            sage: F.prime_above(3, degree=2)
            Traceback (most recent call last):
            ...
            ValueError: No prime of degree 2 above Fractional ideal (3)
            sage: [ id.residue_class_degree() for id, _ in F.ideal(3).factor() ]
            [1]

            Asking for a specific degree works:

            sage: P5_1 = F.prime_above(5, degree=1)
            sage: P5_1 # random
            Fractional ideal (-t^2 - 1)
            sage: P5_1.residue_class_degree()
            1

            sage: P5_2 = F.prime_above(5, degree=2)
            sage: P5_2 # random
            Fractional ideal (t^2 - 2*t - 1)
            sage: P5_2.residue_class_degree()
            2

        TESTS:
            It doesn't make sense to factor the ideal (0):

            sage: F.prime_above(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'NumberFieldIdeal' object has no attribute 'factor'

            Sage can't factor ideals over extension fields yet:

            sage: G = F.extension(x^2 - 11, 'b')
            sage: G.prime_above(13)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        ids = self.primes_above(x, degree)
        if not ids:
            raise ValueError, "No prime of degree %s above %s" % (degree, self.ideal(x))
        return ids[0]

    def primes_of_degree_one_iter(self, num_integer_primes=10000, max_iterations=100):
        r"""
        Return an iterator yielding prime ideals of absolute degree one and small norm.

        WARNING:
            It is possible that there are no primes of $K$ of absolute degree
            one of small prime norm, and it possible that this algorithm will
            not find any primes of small norm.

            See module sage.rings.number_field.small_primes_of_degree_one for details.

        INPUT:
            num_integer_primes (default: 10000) -- an integer.  We try
                to find primes of absolute norm no greater than the
                num_integer_primes-th prime number.  For example, if
                num_integer_primes is 2, the largest norm found will
                be 3, since the second prime is 3.
            max_iterations (default: 100) -- an integer.  We test
                max_iterations integers to find small primes before
                raising StopIteration.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(10)
            sage: it = K.primes_of_degree_one_iter()
            sage: Ps = [ it.next() for i in range(3) ]
            sage: Ps # random
            [Fractional ideal (z^3 + z + 1), Fractional ideal (3*z^3 - z^2 + z - 1), Fractional ideal (2*z^3 - 3*z^2 + z - 2)]
            sage: [ P.norm() for P in Ps ] # random
            [11, 31, 41]
            sage: [ P.residue_class_degree() for P in Ps ]
            [1, 1, 1]
        """
        from sage.rings.number_field.small_primes_of_degree_one import Small_primes_of_degree_one_iter
        return Small_primes_of_degree_one_iter(self, num_integer_primes, max_iterations)

    def primes_of_degree_one_list(self, n, num_integer_primes=10000, max_iterations=100):
        r"""
        Return a list of n prime ideals of absolute degree one and small norm.

        WARNING:
            It is possible that there are no primes of $K$ of absolute degree
            one of small prime norm, and it possible that this algorithm will
            not find any primes of small norm.

            See module sage.rings.number_field.small_primes_of_degree_one for details.

        INPUT:
            num_integer_primes (default: 10000) -- an integer.  We try
                to find primes of absolute norm no greater than the
                num_integer_primes-th prime number.  For example, if
                num_integer_primes is 2, the largest norm found will
                be 3, since the second prime is 3.
            max_iterations (default: 100) -- an integer.  We test
                max_iterations integers to find small primes before
                raising StopIteration.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(10)
            sage: Ps = K.primes_of_degree_one_list(3)
            sage: Ps
            [Fractional ideal (z^3 + z + 1), Fractional ideal (3*z^3 - z^2 + z - 1), Fractional ideal (2*z^3 - 3*z^2 + z - 2)]
            sage: [ P.norm() for P in Ps ]
            [11, 31, 41]
            sage: [ P.residue_class_degree() for P in Ps ]
            [1, 1, 1]
        """
        it = self.primes_of_degree_one_iter()
        return [ it.next() for i in range(n) ]

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        Return whether or not there is a homomorphism defined by the
        given images of generators.

        To do this we just check that the elements of the image of the
        given generator (im_gens always has length 1) satisfies the
        relation of the defining poly of this field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 - 3)
            sage: k._is_valid_homomorphism_(QQ, [0])
            False
            sage: k._is_valid_homomorphism_(k, [])
            False
            sage: k._is_valid_homomorphism_(k, [a])
            True
            sage: k._is_valid_homomorphism_(k, [-a])
            True
            sage: k._is_valid_homomorphism_(k, [a+1])
            False
        """
        try:
            if len(im_gens) != 1:
                return False
            # We need that elements of the base ring of the polynomial
            # ring map canonically into codomain.
            codomain._coerce_(rational.Rational(1))
            f = self.defining_polynomial()
            return codomain(f(im_gens[0])) == 0
        except (TypeError, ValueError):
            return False

    def pari_polynomial(self, name='x'):
        """
        PARI polynomial corresponding to polynomial that defines
        this field. By default, this is a polynomial in the
        variable "x".

        EXAMPLES:
            sage: y = polygen(QQ)
            sage: k.<a> = NumberField(y^2 - 3/2*y + 5/3)
            sage: k.pari_polynomial()
            x^2 - 3/2*x + 5/3
            sage: k.pari_polynomial('a')
            a^2 - 3/2*a + 5/3
        """
        try:
            if (self.__pari_polynomial_var == name):
                return self.__pari_polynomial
            else:
                self.__pari_polynomial = self.__pari_polynomial(name)
                self.__pari_polynomial_var = name
                return self.__pari_polynomial
        except AttributeError:
            poly = self.polynomial()
            with localvars(poly.parent(), name):
                self.__pari_polynomial = poly._pari_()
                self.__pari_polynomial_var = name
            return self.__pari_polynomial

    def pari_nf(self):
        """
        PARI number field corresponding to this field.

        This is the number field constructed using nfinit.
        This is the same as the number field got by doing
        pari(self) or gp(self).

        EXAMPLES:
            sage: k.<a> = NumberField(x^4 - 3*x + 7); k
            Number Field in a with defining polynomial x^4 - 3*x + 7
            sage: k.pari_nf()[:4]
            [x^4 - 3*x + 7, [0, 2], 85621, 1]
            sage: pari(k)[:4]
            [x^4 - 3*x + 7, [0, 2], 85621, 1]

            sage: k.<a> = NumberField(x^4 - 3/2*x + 5/3); k
            Number Field in a with defining polynomial x^4 - 3/2*x + 5/3
            sage: k.pari_nf()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce number field defined by non-integral polynomial to PARI.
            sage: pari(k)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce number field defined by non-integral polynomial to PARI.
            sage: gp(k)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce number field defined by non-integral polynomial to PARI.
        """
        if self.defining_polynomial().denominator() != 1:
            raise TypeError, "Unable to coerce number field defined by non-integral polynomial to PARI."
        try:
            return self.__pari_nf
        except AttributeError:
            f = self.pari_polynomial()
            self.__pari_nf = f.nfinit()
            return self.__pari_nf

    def _pari_init_(self):
        """
        Needed for conversion of number field to PARI.

        This only works if the defining polynomial of this number
        field is integral and monic.

        EXAMPLES:
            sage: k = NumberField(x^2 + x + 1, 'a')
            sage: k._pari_init_()
            'nfinit(x^2 + x + 1)'
            sage: k._pari_()
            [x^2 + x + 1, [0, 1], -3, 1, ... [1, x], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
            sage: pari(k)
            [x^2 + x + 1, [0, 1], -3, 1, ...[1, x], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
        """
        if self.defining_polynomial().denominator() != 1:
            raise TypeError, "Unable to coerce number field defined by non-integral polynomial to PARI."
        return 'nfinit(%s)'%self.pari_polynomial()

    def pari_bnf(self, certify=False, units=True):
        """
        PARI big number field corresponding to this field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: len(k.pari_bnf())
            10
            sage: k.pari_bnf()[:4]
            [[;], matrix(0,7), [;], ...]
            sage: len(k.pari_nf())
            9
        """
        if self.defining_polynomial().denominator() != 1:
            raise TypeError, "Unable to coerce number field defined by non-integral polynomial to PARI."
        try:
            if certify:
                self.pari_bnf_certify()
            return self.__pari_bnf
        except AttributeError:
            f = self.pari_polynomial()
            if units:
                self.__pari_bnf = f.bnfinit(1)
            else:
                self.__pari_bnf = f.bnfinit()
            if certify:
                self.pari_bnf_certify()
            return self.__pari_bnf

    def pari_bnf_certify(self):
        """
        Run the PARI bnfcertify function to ensure the correctness of answers.

        If this function returns True (and doesn't raise a
        ValueError), then certification succeeded, and results that
        use the PARI bnf structure with this field are supposed to be
        correct.

        WARNING: I wouldn't trust this to mean that everything
        computed involving this number field is actually correct.

        EXAMPLES:
            sage: k.<a> = NumberField(x^7 + 7); k
            Number Field in a with defining polynomial x^7 + 7
            sage: k.pari_bnf_certify()
            True
        """
        if self.defining_polynomial().denominator() != 1:
            raise TypeError, "Unable to coerce number field defined by non-integral polynomial to PARI."
        if not self.__pari_bnf_certified:
            if self.pari_bnf(certify=False, units=True).bnfcertify() != 1:
                raise ValueError, "The result is not correct according to bnfcertify"
            self.__pari_bnf_certified = True
        return self.__pari_bnf_certified

    def _gap_init_(self):
        """
        Create a gap object representing self and return its name

        EXAMPLE:
            sage: F=CyclotomicField(8)
            sage: F.gen()
            zeta8
            sage: F._gap_init_() # the following variable name $sage1 represents the F.base_ring() in gap and is somehow random
            'CallFuncList(function() local x,E; x:=Indeterminate($sage1,"x"); E:=AlgebraicExtension($sage1,x^4 + 1,"zeta8"); return E; end,[])'
            sage: f=gap(F)
            sage: f
            <algebraic extension over the Rationals of degree 4>
            sage: f.GeneratorsOfDivisionRing()
            [ (zeta8) ]
        """
        if not self.is_absolute():
            raise NotImplementedError, "Currently, only simple algebraic extensions are implemented in gap"
        G = sage.interfaces.gap.gap
        return 'CallFuncList(function() local x,E; x:=Indeterminate(%s,"x"); E:=AlgebraicExtension(%s,%s,"%s"); return E; end,[])'%(G(self.base_ring()).name(),G(self.base_ring()).name(),self.polynomial().__repr__(),str(self.gen()))

    def characteristic(self):
        """
        Return the characteristic of this number field, which is
        of course 0.

        EXAMPLES:
            sage: k.<a> = NumberField(x^99 + 2); k
            Number Field in a with defining polynomial x^99 + 2
            sage: k.characteristic()
            0
        """
        return 0

    def class_group(self, proof=None, names='c'):
        r"""
        Return the class group of the ring of integers of this number field.

        INPUT:
            proof -- if True then compute the classgroup provably correctly.
                     Default is True.  Call number_field_proof to change
                     this default globally.
            names -- names of the generators of this class group.

        OUTPUT:
            The class group of this number field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 23)
            sage: G = K.class_group(); G
            Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
            sage: G.0
            Fractional ideal class (2, 1/2*a - 1/2)
            sage: G.gens()
            [Fractional ideal class (2, 1/2*a - 1/2)]

            sage: G.number_field()
            Number Field in a with defining polynomial x^2 + 23
            sage: G is K.class_group()
            True
            sage: G is K.class_group(proof=False)
            False
            sage: G.gens()
            [Fractional ideal class (2, 1/2*a - 1/2)]

        There can be multiple generators:
            sage: k.<a> = NumberField(x^2 + 20072)
            sage: G = k.class_group(); G
            Class group of order 76 with structure C38 x C2 of Number Field in a with defining polynomial x^2 + 20072
            sage: G.gens()
            [Fractional ideal class (41, a + 10), Fractional ideal class (2, -1/2*a)]
            sage: G.0
            Fractional ideal class (41, a + 10)
            sage: G.0^20
            Fractional ideal class (43, a + 3)
            sage: G.0^38
            Trivial principal fractional ideal class
            sage: G.1
            Fractional ideal class (2, -1/2*a)
            sage: G.1^2
            Trivial principal fractional ideal class

        Class groups of Hecke polynomials tend to be very small:
            sage: f = ModularForms(97, 2).T(2).charpoly()
            sage: f.factor()
            (x - 3) * (x^3 + 4*x^2 + 3*x - 1) * (x^4 - 3*x^3 - x^2 + 6*x - 1)
            sage: for g,_ in f.factor(): print NumberField(g,'a').class_group().order()
            ...
            1
            1
            1
        """
        proof = proof_flag(proof)
        try:
            return self.__class_group[proof, names]
        except KeyError:
            pass
        except AttributeError:
            self.__class_group = {}
        k = self.pari_bnf(proof)
        cycle_structure = eval(str(k.getattr('clgp.cyc')))

        # First gens is a pari list of pari gens
        gens = k.getattr('clgp.gen')
        R    = self.polynomial_ring()

        # Next gens is a list of ideals.
        gens = [self.ideal([self(R(convert_from_zk_basis(self, y))) for y in x]) for x in gens]

        G    = ClassGroup(cycle_structure, names, self, gens)
        self.__class_group[proof,names] = G
        return G

    def class_number(self, proof=None):
        """
        Return the class number of this number field, as an integer.

        INPUT:
            proof -- bool (default: True unless you called number_field_proof)

        EXAMPLES:
            sage: NumberField(x^2 + 23, 'a').class_number()
            3
            sage: NumberField(x^2 + 163, 'a').class_number()
            1
            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').class_number(proof=False)
            1539
        """
        proof = proof_flag(proof)
        return self.class_group(proof).order()

    def composite_fields(self, other, names=None):
        """
        List of all possible composite number fields formed from self
        and other.

        INPUT:
            other -- a number field
            names -- generator name for composite fields

        OUTPUT:
            list -- list of the composite fields.

        EXAMPLES:
            sage: K.<a> = NumberField(x^4 - 2)
            sage: K.composite_fields(K)
            [Number Field in a0 with defining polynomial x^4 - 162,
             Number Field in a1 with defining polynomial x^4 - 2,
             Number Field in a2 with defining polynomial x^8 + 28*x^4 + 2500]
            sage: k.<a> = NumberField(x^3 + 2)
            sage: m.<b> = NumberField(x^3 + 2)
            sage: k.composite_fields(m, 'c')
            [Number Field in c0 with defining polynomial x^3 - 2,
             Number Field in c1 with defining polynomial x^6 - 40*x^3 + 1372]
        """
        if names is None:
            sv = self.variable_name(); ov = other.variable_name()
            names = sv + (ov if ov != sv else "")
        if not isinstance(other, NumberField_generic):
            raise TypeError, "other must be a number field."
        f = self.pari_polynomial()
        g = other.pari_polynomial()
        C = f.polcompositum(g)
        R = self.polynomial().parent()
        C = [R(h) for h in C]
        if len(C) == 1:
            return [NumberField(C[0], names)]
        else:
            name = sage.structure.parent_gens.normalize_names(1, names)[0]
            return [NumberField(C[i], name + str(i)) for i in range(len(C))]

    def absolute_degree(self):
        """
        Return the degree of self over $\QQ$.

        EXAMPLES:
            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').absolute_degree()
            3
            sage: NumberField(x + 1, 'a').absolute_degree()
            1
            sage: NumberField(x^997 + 17*x + 3, 'a', check=False).absolute_degree()
            997
        """
        return self.polynomial().degree()

    def degree(self):
        """
        Return the degree of this number field.

        EXAMPLES:
            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').degree()
            3
            sage: NumberField(x + 1, 'a').degree()
            1
            sage: NumberField(x^997 + 17*x + 3, 'a', check=False).degree()
            997
        """
        return self.polynomial().degree()

    def different(self):
        r"""
        Compute the different fractional ideal of this number field.

        The different is the set of all $x$ in $K$ such that the trace
        of $xy$ is an integer for all $y \in O_K$.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 23)
            sage: d = k.different()
            sage: d        # random sign in output
            Fractional ideal (-a)
            sage: d.norm()
            23
            sage: k.disc()
            -23

        The different is cached:
            sage: d is k.different()
            True

        Another example:
            sage: k.<b> = NumberField(x^2 - 123)
            sage: d = k.different(); d
            Fractional ideal (2*b)
            sage: d.norm()
            492
            sage: k.disc()
            492
        """
        try:
            return self.__different

        except AttributeError:

            diff = self.pari_nf().getattr('diff')
            zk_basis = self.pari_nf().getattr('zk')
            basis_elts = zk_basis * diff
            R = self.polynomial().parent()
            self.__different = self.ideal([ self(R(x)) for x in basis_elts ])
            return self.__different

    def discriminant(self, v=None):
        """
        Returns the discriminant of the ring of integers of the number field,
        or if v is specified, the determinant of the trace pairing
        on the elements of the list v.

        INPUT:
            v (optional) -- list of element of this number field
        OUTPUT:
            Integer if v is omitted, and Rational otherwise.

        EXAMPLES:
            sage: K.<t> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: K.disc()
            -503
            sage: K.disc([1, t, t^2])
            -2012
            sage: K.disc([1/7, (1/5)*t, (1/3)*t^2])
            -2012/11025
            sage: (5*7*3)^2
            11025
        """
        if v == None:
            try:
                return self.__disc
            except AttributeError:
                self.__disc = ZZ(str(self.pari_polynomial().nfdisc()))
                return self.__disc
        else:
            return QQ(self.trace_pairing(v).det())

    def disc(self, v=None):
        """
        Shortcut for self.discriminant.

        EXAMPLES:
            sage: k.<b> = NumberField(x^2 - 123)
            sage: k.disc()
            492
        """
        return self.discriminant(v=v)

    def elements_of_norm(self, n, proof=None):
        r"""
        Return a list of solutions modulo units of positive norm to
        $Norm(a) = n$, where a can be any integer in this number field.

        INPUT:
            proof -- default: True, unless you called number_field_proof and
            set it otherwise.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2+1)
            sage: K.elements_of_norm(3)
            []
            sage: K.elements_of_norm(50)
            [7*a - 1, -5*a + 5, a - 7]           # 32-bit
            [7*a - 1, -5*a + 5, -7*a - 1]        # 64-bit
        """
        proof = proof_flag(proof)
        B = self.pari_bnf(proof).bnfisintnorm(n)
        R = self.polynomial().parent()
        return [self(QQ['x'](R(g))) for g in B]

    def extension(self, poly, name=None, names=None, check=True):
        """
        Return the relative extension of this field by a given polynomial.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 - 2)
            sage: R.<t> = K[]
            sage: L.<b> = K.extension(t^2 + a); L
            Number Field in b with defining polynomial t^2 + a over its base field

        We create another extension.
            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: y = var('y')
            sage: m.<b> = k.extension(y^2 + 2); m
            Number Field in b with defining polynomial y^2 + 2 over its base field

        Note that b is a root of $y^2 + 2$:
            sage: b.minpoly()
            x^2 + 2
            sage: b.minpoly('z')
            z^2 + 2

        A relative extension of a relative extension.
            sage: k.<a> = NumberField([x^2 + 1, x^3 + x + 1])
            sage: R.<z> = k[]
            sage: L.<b> = NumberField(z^3 + 3 + a); L
            Number Field in b with defining polynomial z^3 + a0 + 3 over its base field
        """
        if not isinstance(poly, polynomial_element.Polynomial):
            try:
                poly = poly.polynomial(self)
            except (AttributeError, TypeError):
                raise TypeError, "polynomial (=%s) must be a polynomial."%repr(poly)
        if not names is None:
            name = names
        if isinstance(name, tuple):
            name = name[0]
        if name is None:
            raise TypeError, "the variable name must be specified."
        return NumberField_relative(self, poly, str(name), check=check)

    def factor(self, n):
        r"""
        Ideal factorization of the principal ideal of the ring
        of integers generated by $n$.

        EXAMPLE:
        Here we show how to factor gaussian integers.
        First we form a number field defined by $x^2 + 1$:

            sage: K.<I> = NumberField(x^2 + 1); K
            Number Field in I with defining polynomial x^2 + 1

        Here are the factors:

            sage: fi, fj = K.factor(13); fi,fj
            ((Fractional ideal (-3*I - 2), 1), (Fractional ideal (3*I - 2), 1))

        Now we extract the reduced form of the generators:

            sage: zi = fi[0].gens_reduced()[0]; zi
            -3*I - 2
            sage: zj = fj[0].gens_reduced()[0]; zj
            3*I - 2

        We recover the integer that was factored in $\Z[i]$

            sage: zi*zj
            13

        One can also factor elements of the number field:

            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.factor(1/3)
            Fractional ideal (3)^-1
            sage: K.factor(1+a)
            Fractional ideal (a + 1)
            sage: K.factor(1+a/5)
            (Fractional ideal (-3*a - 2)) * (Fractional ideal (a + 1)) * (Fractional ideal (-a - 2))^-1 * (Fractional ideal (2*a + 1))^-1

        AUTHOR:
            -- Alex Clemesha (2006-05-20): examples
        """
        return self.ideal(n).factor()

    def gen(self, n=0):
        """
        Return the generator for this number field.

        INPUT:
            n -- must be 0 (the default), or an exception is raised.

        EXAMPLES:
            sage: k.<theta> = NumberField(x^14 + 2); k
            Number Field in theta with defining polynomial x^14 + 2
            sage: k.gen()
            theta
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: Only one generator.
        """
        if n != 0:
            raise IndexError, "Only one generator."
        try:
            return self.__gen
        except AttributeError:
            if self.__polynomial != None:
                X = self.__polynomial.parent().gen()
            else:
                X = PolynomialRing(rational_field.RationalField()).gen()
            self.__gen = self._element_class(self, X)
            return self.__gen

    def is_field(self):
        """
        Return True since a number field is a field.

        EXAMPLES:
            sage: NumberField(x^5 + x + 3, 'c').is_field()
            True
        """
        return True

    def is_galois(self):
        r"""
        Return True if this number field is a Galois extension of $\QQ$.

        EXAMPLES:
            sage: NumberField(x^2 + 1, 'i').is_galois()
            True
            sage: NumberField(x^3 + 2, 'a').is_galois()
            False
        """
        return self.galois_group(pari_group=True).group().order() == self.degree()

    def galois_group(self, pari_group = True, algorithm='pari'):
        r"""
        Return the Galois group of the Galois closure of this number
        field as an abstract group.

        INPUT:
            pari_group -- bool (default: False); if True instead return
                          the Galois group as a PARI group.
            algorithm -- 'pari', 'kash', 'magma' (default: 'pari', except
                          when the degree is >= 12 when 'kash' is tried)

        For more (important!) documentation, so the documentation
        for Galois groups of polynomials over $\QQ$, e.g., by
        typing \code{K.polynomial().galois_group?}, where $K$
        is a number field.

        To obtain actual field automorphisms that can be applied to
        elements, use \code{End(K).list()} and
        \code{K.galois_closure()} together (see example below).

        EXAMPLES:
            sage: k.<b> = NumberField(x^2 - 14)
            sage: k.galois_group ()
            Galois group PARI group [2, -1, 1, "S2"] of degree 2 of the Number Field in b with defining polynomial x^2 - 14

            sage: NumberField(x^3-2, 'a').galois_group(pari_group=True)
            Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in a with defining polynomial x^3 - 2

            sage: NumberField(x-1, 'a').galois_group(pari_group=False)    # optional - database_gap
            Galois group Transitive group number 1 of degree 1 of the Number Field in a with defining polynomial x - 1
            sage: NumberField(x^2+2, 'a').galois_group(pari_group=False)  # optional - database_gap
            Galois group Transitive group number 1 of degree 2 of the Number Field in a with defining polynomial x^2 + 2
            sage: NumberField(x^3-2, 'a').galois_group(pari_group=False)  # optional - database_gap
            Galois group Transitive group number 2 of degree 3 of the Number Field in a with defining polynomial x^3 - 2

            sage: x = polygen(QQ)
            sage: NumberField(x^3 + 2*x + 1, 'a').galois_group(pari_group=False)
            Galois group Transitive group number 2 of degree 3 of the Number Field in a with defining polynomial x^3 + 2*x + 1
            sage: NumberField(x^3 + 2*x + 1, 'a').galois_group(algorithm='magma')   # optional - magma
            Galois group Transitive group number 2 of degree 3 of the Number Field in a with defining polynomial x^3 + 2*x + 1



        EXPLICIT GALOIS GROUP:
        We compute the Galois group as an explicit group of
        automorphisms of the Galois closure of a field.

            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<b1> = K.galois_closure(); L
            Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: G = End(L); G
            Automorphism group of Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: G.list()
            [
            Ring endomorphism of Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
              Defn: b1 |--> b1,
            ...
            Ring endomorphism of Number Field in b1 with defining polynomial x^6 + 40*x^3 + 1372
              Defn: b1 |--> -2/63*b1^4 - 31/63*b1
            ]
            sage: G[1](b1)
            1/36*b1^4 + 1/18*b1
        """
        try:
            return self.__galois_group[pari_group, algorithm]
        except KeyError:
            pass
        except AttributeError:
            self.__galois_group = {}

        G = self.polynomial().galois_group(pari_group = pari_group, algorithm = algorithm)
        H = GaloisGroup(G, self)
        self.__galois_group[pari_group, algorithm] = H
        return H

    def _normalize_prime_list(self, v):
        """
        Internal function to convert into a tuple of primes either None
        or a single prime or a list.

        EXAMPLES:
            sage: K.<i> = NumberField(x^2 + 1)
            sage: K._normalize_prime_list(None)
            ()
            sage: K._normalize_prime_list(3)
            (3,)
            sage: K._normalize_prime_list([3,5])
            (3, 5)
        """
        if v is None:
            v = []
        elif not isinstance(v, (list, tuple)):
            v = [v]
        return tuple([ZZ(x) for x in v])

    def power_basis(self):
        r"""
        Return a power basis for this number field over its base field.

        If this number field is represented as $k[t]/f(t)$, then the basis
        returned is $1, t, t^2, \ldots, t^{d-1}$ where $d$ is the degree of
        this number field over its base field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 + 10*x + 1)
            sage: K.power_basis()
            [1, a, a^2, a^3, a^4]

            sage: L.<b> = K.extension(x^2 - 2)
            sage: L.power_basis()
            [1, b]
            sage: L.absolute_field('c').power_basis()
            [1, c, c^2, c^3, c^4, c^5, c^6, c^7, c^8, c^9]

            sage: M = CyclotomicField(15)
            sage: M.power_basis()
            [1, zeta15, zeta15^2, zeta15^3, zeta15^4, zeta15^5, zeta15^6, zeta15^7]
        """
        g = self.gen()
        return [ g**i for i in range(self.degree()) ]

    def integral_basis(self, v=None):
        """
        Returns a list containing a ZZ-basis for the full ring of
        integers of this number field.

        INPUT:
            v -- None, a prime, or a list of primes.  See
            the documentation for self.maximal_order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 + 10*x + 1)
            sage: K.integral_basis()
            [1, a, a^2, a^3, a^4]

        Next we compute the ring of integers of a cubic field in which 2
        is an "essential discriminant divisor", so the ring of integers
        is not generated by a single element.
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: K.integral_basis()
            [1, 1/2*a^2 + 1/2*a, a^2]

        ALGORITHM:
            Uses the pari library.
        """
        return self.maximal_order().basis()

    def _compute_integral_basis(self, v=None):
        """
        Internal function returning an integral basis of this number
        field; used in the maximal_order() function.

        Note that this is *not* necessarily the same basis
        returned by self.integral_basis().

        INPUT:
            v -- None, a prime, or a list of primes.  See
            the documentation for self.maximal_order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 + 10*x + 1)
            sage: K._compute_integral_basis()
            [1, a, a^2, a^3, a^4]

        Next we compute the ring of integers of a cubic field in which 2
        is an "essential discriminant divisor", so the ring of integers
        is not generated by a single element.
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: K._compute_integral_basis()
            [1, a, 1/2*a^2 + 1/2*a]
            sage: K.integral_basis()
            [1, 1/2*a^2 + 1/2*a, a^2]
        """
        v = self._normalize_prime_list(v)
        try:
            return self._integral_basis_dict[v]
        except:
            f = self.pari_polynomial()

            if len(v) == 0:
                B = f.nfbasis()
            else:
                m = pari.matrix(len(v), 2)
                d = f.poldisc()
                for i in range(len(v)):
                    p = pari(ZZ(v[i]))
                    m[i,0] = p
                    m[i,1] = d.valuation(p)
                B = f.nfbasis(p = m)

            R = self.polynomial().parent()
            basis = [self(R(g).list()) for g in B]
            self._integral_basis_dict[v] = basis
            return basis

    def reduced_basis(self, prec=None):
        r"""
        This function returns an LLL-reduced basis for the
        Minkowski-embedding of the maximal order of a number field.

        INPUT:
            self -- number field, the base field
            prec (default: None) -- the precision with which
              to compute the Minkowski embedding. (See NOTE below.)

        OUTPUT:
            An LLL-reduced basis for the Minkowski-embedding of the
            maximal order of a number field, given by a sequence of
            (integral) elements from the field.

        NOTE: In the non-totally-real case, the LLL routine we call is
            currently Pari's qflll(), which works with floating point
            approximations, and so the result is only as good as the
            precision promised by Pari. The matrix returned will
            always be integral; however, it may only be only "almost"
            LLL-reduced when the precision is not sufficiently high.

            If the following run-time error occurs:
            "PariError: not a definite matrix in lllgram (42)"
            try increasing the prec parameter,

        EXAMPLES:
            sage: F.<t> = NumberField(x^6-7*x^4-x^3+11*x^2+x-1)
            sage: F.maximal_order().basis()
            [1/2*t^5 + 1/2*t^4 + 1/2*t^2 + 1/2, t, t^2, t^3, t^4, t^5]
            sage: F.reduced_basis()
            [-1, -1/2*t^5 + 1/2*t^4 + 3*t^3 - 3/2*t^2 - 4*t - 1/2, t, 1/2*t^5 + 1/2*t^4 - 4*t^3 - 5/2*t^2 + 7*t + 1/2, 1/2*t^5 - 1/2*t^4 - 2*t^3 + 3/2*t^2 - 1/2, 1/2*t^5 - 1/2*t^4 - 3*t^3 + 5/2*t^2 + 4*t - 5/2]

            sage: F.<alpha> = NumberField(x^4+x^2+712312*x+131001238)
            sage: F.integral_basis()
            [1, alpha, 1/2*alpha^3 + 1/2*alpha^2, alpha^3]
            sage: F.reduced_basis(prec=64)
            [1, alpha, alpha^3 - 2*alpha^2 + 15*alpha, 16*alpha^3 - 31*alpha^2 + 469*alpha + 267109]     # 32-bit
            Traceback (most recent call last):                 # 64-bit
            ...                                                # 64-bit
            PariError: not a definite matrix in lllgram (42)   # 64-bit
            sage: F.reduced_basis(prec=96)
            [1, alpha, alpha^3 - 2*alpha^2 + 15*alpha, 16*alpha^3 - 31*alpha^2 + 469*alpha + 267109]
        """
        if self.is_totally_real():
            try:
                return self.__reduced_basis
            except AttributeError:
                pass
        else:
            try:
                if self.__reduced_basis_precision >= prec:
                    return self.__reduced_basis
            except AttributeError:
                pass

        from sage.matrix.constructor import matrix

        d = self.degree()
        Z_basis = self.integral_basis()

        ## If self is totally real, then we can use (x*y).trace() as
        ## the inner product on the Minkowski embedding, which is
        ## faster than computing all the conjugates, etc ...
        if self.is_totally_real():
            T = pari(matrix(ZZ, d, d, [[(x*y).trace() for x in Z_basis]
                                       for y in Z_basis])).qflllgram()
            self.__reduced_basis = [ sum([ ZZ(T[i][j]) * Z_basis[j]
                                           for j in range(d)])
                                     for i in range(d)]
        else:
            M = self.Minkowski_embedding(self.integral_basis(), prec=prec)
            T = pari(M).qflll().python()
            self.__reduced_basis = [ self(v) for v in T.columns() ]
            if prec is None:
                ## this is the default choice for Minkowski_embedding
                self.__reduced_basis_prec = 53
            else:
                self.__reduced_basis_prec = prec

        return self.__reduced_basis


    def reduced_gram_matrix(self, prec=None):
        r"""
        This function returns the Gram matrix of an LLL-reduced basis
        for the Minkowski embedding of the maximal order of a number
        field.

        INPUT:
            self -- number field, the base field
            prec (default: None) -- the precision with which to
              calculate the Minkowski embedding. (See NOTE below.)

        OUTPUT:
            The Gram matrix $[<x_i,x_j>]$ of an LLL reduced basis for
            the maximal order of self, where the integral basis for
            self is given by $\{x_0, \dots, x_{n-1}\}$. Here < , > is
            the usual inner product on $\RR^n$, and self is
            embedded in $\RR^n$ by the Minkowski embedding. See
            the docstring for self.Minkowski_embedding for more
            information.

        NOTE: In the non-totally-real case, the LLL routine we call is
            currently Pari's qflll(), which works with floating point
            approximations, and so the result is only as good as the
            precision promised by Pari.  In particular, in this case,
            the returned matrix will *not* be integral, and may not
            have enough precision to recover the correct gram matrix
            (which is known to be integral for theoretical reasons).
            Thus the need for the prec flag above.

            If the following run-time error occurs:
            "PariError: not a definite matrix in lllgram (42)"
            try increasing the prec parameter,

        EXAMPLES:
            sage: F.<t> = NumberField(x^6-7*x^4-x^3+11*x^2+x-1)
            sage: F.reduced_gram_matrix()
            [ 6  3  0  2  0  1]
            [ 3  9  0  1  0 -2]
            [ 0  0 14  6 -2  3]
            [ 2  1  6 16 -3  3]
            [ 0  0 -2 -3 16  6]
            [ 1 -2  3  3  6 19]
            sage: Matrix(6, [(x*y).trace() for x in F.integral_basis() for y in F.integral_basis()])
            [2550  133  259  664 1368 3421]
            [ 133   14    3   54   30  233]
            [ 259    3   54   30  233  217]
            [ 664   54   30  233  217 1078]
            [1368   30  233  217 1078 1371]
            [3421  233  217 1078 1371 5224]

            sage: var('x')
            x
            sage: F.<alpha> = NumberField(x^4+x^2+712312*x+131001238)
            sage: F.reduced_gram_matrix(prec=128)
            [   4.0000000000000000000000000000000000000   0.00000000000000000000000000000000000000 -2.1369320000000000000000000000000000000e6 -3.3122478000000000000000000000000000000e7]
            [  0.00000000000000000000000000000000000000    46721.539331563218381658483353092335550 -2.2467769057394530109094755223395819322e7 -3.4807276041138450473611629088647496430e8]
            [-2.1369320000000000000000000000000000000e6 -2.2467769057394530109094755223395819322e7 7.0704243186034907491782135494859351061e12 1.1256636615786237006027526953641297995e14]
            [-3.3122478000000000000000000000000000000e7 -3.4807276041138450473611629088647496430e8 1.1256636615786237006027526953641297995e14 1.7923838231014970520503146603069479547e15]
        """
        if self.is_totally_real():
            try:
                return self.__reduced_gram_matrix
            except AttributeError:
                pass
        else:
            try:
                if self.__reduced_gram_matrix_prec >= prec:
                    return self.__reduced_gram_matrix
            except AttributeError:
                pass

        from sage.matrix.constructor import matrix
        from sage.misc.flatten import flatten
        d = self.degree()

        if self.is_totally_real():
            B = self.reduced_basis()
            self.__reduced_gram_matrix = matrix(ZZ, d, d,
                                                [[(x*y).trace() for x in B]
                                                 for y in B])
        else:
            M = self.Minkowski_embedding(prec=prec)
            T = matrix(d, flatten([ a.vector().list()
                                    for a in self.reduced_basis(prec=prec) ]))
            A = M*(T.transpose())
            self.__reduced_gram_matrix = A.transpose()*A
            if prec is None:
                ## this is the default choice for Minkowski_embedding
                self.__reduced_gram_matrix_prec = 53
            else:
                self.__reduced_gram_matrix_prec = prec

        return self.__reduced_gram_matrix


    #******************************************************
    # Supplementary algorithm to enumerate lattice points
    #******************************************************

    def _positive_integral_elements_with_trace(self, C):

        r"""
        Find all totally positive integral elements in self whose
        trace is between C[0] and C[1], inclusive.

        NOTE: This is currently only implemented in the case
        that self is totally real, since it requires exact
        computation of self.reduced_gram_matrix().

        EXAMPLES:
            sage: K.<alpha> = NumberField(ZZ['x'].0^2-2)
            sage: K._positive_integral_elements_with_trace([0,5])
            [alpha + 2, -alpha + 2, 2, 1]
            sage: L.<beta> = NumberField(ZZ['x'].0^2+1)
            sage: L._positive_integral_elements_with_trace([5,11])
            Traceback (most recent call last):
            ...
            NotImplementedError: exact computation of LLL reduction only implemented in the totally real case
            sage: L._positive_integral_elements_with_trace([-5,1])
            Traceback (most recent call last):
            ...
            ValueError: bounds must be positive
        """
        if C[0] < 0:
            raise ValueError, "bounds must be positive"

        if not self.is_totally_real():
            raise NotImplementedError, "exact computation of LLL reduction only implemented in the totally real case"

        Z_F = self.maximal_order()
        B = self.reduced_basis()
        T = self.reduced_gram_matrix()
        P = pari(T).qfminim((C[1]**2)*(1./2), 10**6)[2]

        S = []
        for p in P:
            theta = sum([ p.list()[i]*B[i] for i in range(self.degree())])
            if theta.trace() < 0:
                theta *= -1
            if theta.trace() >= C[0] and theta.trace() <= C[1]:
                inbounds = True
                for v in self.real_embeddings():
                    inbounds = inbounds and v(theta) > 0
                if inbounds:
                    S.append(self(theta))
        return S


    def zeta_function(self, prec=53,
                      max_imaginary_part=0,
                      max_asymp_coeffs=40):
        r"""
        Return the Zeta function of this number field.

        This actually returns an interface to Tim Dokchitser's program
        for computing with the Dedekind zeta function zeta_F(s) of the
        number field F.

        INPUT:
            prec -- integer (bits precision)
            max_imaginary_part -- real number
            max_asymp_coeffs -- integer

        OUTPUT:
            The zeta function of this number field.

        EXAMPLES:
            sage: K.<a> = NumberField(ZZ['x'].0^2+ZZ['x'].0-1)
            sage: Z = K.zeta_function()
            sage: Z
            Zeta function associated to Number Field in a with defining polynomial x^2 + x - 1
            sage: Z(-1)
            0.0333333333333333
        """
        from sage.lfunctions.all import Dokchitser
        key = (prec, max_imaginary_part, max_asymp_coeffs)
        r1 = self.signature()[0]
        r2 = self.signature()[1]
        zero = [0]
        one = [1]
        Z = Dokchitser(conductor = abs(self.discriminant()),
                       gammaV = (r1+r2)*zero + r2*one,
                       weight = 1,
                       eps = 1,
                       poles = [1],
                       prec = prec)
        s = 'nf = nfinit(%s);'%self.polynomial()
        s += 'dzk = dirzetak(nf,cflength());'
        Z.init_coeffs('dzk[k]',pari_precode = s,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        Z.check_functional_equation()
        Z.rename('Zeta function associated to %s'%self)
        return Z

    def narrow_class_group(self, proof=None):
        r"""
        Return the narrow class group of this field.

        INPUT:
            proof -- default: None (use the global proof setting, which defaults to True).

        EXAMPLES:
            sage: NumberField(x^3+x+9, 'a').narrow_class_group()
            Multiplicative Abelian Group isomorphic to C2
        """
        proof = proof_flag(proof)
        try:
            return self.__narrow_class_group
        except AttributeError:
            k = self.pari_bnf(proof)
            s = str(k.bnfnarrow())
            s = s.replace(";",",")
            s = eval(s)
            self.__narrow_class_group = sage.groups.abelian_gps.abelian_group.AbelianGroup(s[1])
        return self.__narrow_class_group

    def ngens(self):
        """
        Return the number of generators of this number field (always 1).

        OUTPUT:
            the python integer 1.

        EXAMPLES:
            sage: NumberField(x^2 + 17,'a').ngens()
            1
            sage: NumberField(x + 3,'a').ngens()
            1
            sage: k.<a> = NumberField(x + 3)
            sage: k.ngens()
            1
            sage: k.0
            -3
        """
        return 1

    def order(self):
        """
        Return the order of this number field (always +infinity).

        OUTPUT:
            always positive infinity

        EXAMPLES:
            sage: NumberField(x^2 + 19,'a').order()
            +Infinity
        """
        return infinity.infinity

    def polynomial_ntl(self):
        """
        Return defining polynomial of this number field
        as a pair, an ntl polynomial and a denominator.

        This is used mainly to implement some internal arithmetic.

        EXAMPLES:
            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').polynomial_ntl()
            ([-27 34 51], 51)
        """
        try:
            return (self.__polynomial_ntl, self.__denominator_ntl)
        except AttributeError:
            self.__denominator_ntl = ntl.ZZ()
            den = self.polynomial().denominator()
            self.__denominator_ntl.set_from_sage_int(ZZ(den))
            self.__polynomial_ntl = ntl.ZZX((self.polynomial()*den).list())
        return (self.__polynomial_ntl, self.__denominator_ntl)

    def polynomial(self):
        """
        Return the defining polynomial of this number field.

        This is exactly the same as \code{self.defining_polynomal()}.

        EXAMPLES:
            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').polynomial()
            x^2 + 2/3*x - 9/17
        """
        return self.__polynomial

    def defining_polynomial(self):   # do not overload this -- overload polynomial instead
        """
        Return the defining polynomial of this number field.

        This is exactly the same as \code{self.polynomal()}.

        EXAMPLES:
            sage: k5.<z> = CyclotomicField(5)
            sage: k5.defining_polynomial()
            x^4 + x^3 + x^2 + x + 1
            sage: y = polygen(QQ,'y')
            sage: k.<a> = NumberField(y^9 - 3*y + 5); k
            Number Field in a with defining polynomial y^9 - 3*y + 5
            sage: k.defining_polynomial()
            y^9 - 3*y + 5
        """
        return self.polynomial()

    def polynomial_ring(self):
        """
        Return the polynomial ring that we view this number field as being
        a quotient of (by a principal ideal).

        EXAMPLES:
        An example with an absolute field:
            sage: k.<a> = NumberField(x^2 + 3)
            sage: y = polygen(QQ, 'y')
            sage: k.<a> = NumberField(y^2 + 3)
            sage: k.polynomial_ring()
            Univariate Polynomial Ring in y over Rational Field

        An example with a relative field:
            sage: y = polygen(QQ, 'y')
            sage: M.<a> = NumberField([y^3 + 97, y^2 + 1]); M
            Number Field in a0 with defining polynomial y^3 + 97 over its base field
            sage: M.polynomial_ring()
            Univariate Polynomial Ring in y over Number Field in a1 with defining polynomial y^2 + 1
        """
        return self.polynomial().parent()

    def polynomial_quotient_ring(self):
        """
        Return the polynomial quotient ring isomorphic to this number field.

        EXAMPLES:
            sage: K = NumberField(x^3 + 2*x - 5, 'alpha')
            sage: K.polynomial_quotient_ring()
            Univariate Quotient Polynomial Ring in alpha over Rational Field with modulus x^3 + 2*x - 5
        """
        return self.polynomial_ring().quotient(self.polynomial(), self.variable_name())

    def regulator(self, proof=None):
        """
        Return the regulator of this number field.

        Note that PARI computes the regulator to higher precision than
        the SAGE default.

        INPUT:
            proof -- default: True, unless you set it otherwise.

        EXAMPLES:
            sage: NumberField(x^2-2, 'a').regulator()
            0.881373587019543
            sage: NumberField(x^4+x^3+x^2+x+1, 'a').regulator()
            0.962423650119207
        """
        proof = proof_flag(proof)
        try:
            return self.__regulator
        except AttributeError:
            from sage.rings.all import RealField
            R = RealField(53)
            k = self.pari_bnf(proof)
            s = str(k.getattr('reg'))
            self.__regulator = R(s)
        return self.__regulator

    def residue_field(self, prime, names = None, check = True):
        """
        Return the residue field of this number field at a given prime, ie $O_K / p O_K$.

        INPUT:
            prime -- a prime ideal of the maximal order in this number
                     field, or an element of the field which generates
                     a principal prime ideal.
            names -- the name of the variable in the residue field
            check -- whether or not to check the primality of prime.
        OUTPUT:
            The residue field at this prime.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: K.residue_field(P)
            Residue field in abar of Fractional ideal (-2*a^2 + 1)

            sage: K.<i> = NumberField(x^2 + 1)
            sage: K.residue_field(1+i)
            Residue field of Fractional ideal (i + 1)

        TESTS:
            sage: L.<b> = NumberField(x^2 + 5)
            sage: L.residue_field(P)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (-2*a^2 + 1) is not an ideal of Number Field in b with defining polynomial x^2 + 5
            sage: L.residue_field(2)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (2) is not a prime ideal

            sage: L.residue_field(L.prime_above(5)^2)
            Traceback (most recent call last):
            ...
            ValueError: Fractional ideal (5) is not a prime ideal
        """
        # This allows principal ideals to be specified using a generator:
        try:
            prime = self.ideal(prime)
        except TypeError:
            pass

        from sage.rings.number_field.number_field_ideal import is_NumberFieldIdeal
        if not is_NumberFieldIdeal(prime) or prime.number_field() is not self:
            raise ValueError, "%s is not an ideal of %s"%(prime,self)
        if check and not prime.is_prime():
            raise ValueError, "%s is not a prime ideal"%prime
        from sage.rings.residue_field import ResidueField
        return ResidueField(prime, names = names, check = False)

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.

        EXAMPLES:
            sage: NumberField(x^2+1, 'a').signature()
            (0, 1)
            sage: NumberField(x^3-2, 'a').signature()
            (1, 1)
        """
        r1, r2 = self.pari_nf().getattr('sign')
        return (ZZ(r1), ZZ(r2))

    def trace_pairing(self, v):
        """
        Return the matrix of the trace pairing on the elements of the
        list $v$.

        EXAMPLES:
            sage: K.<zeta3> = NumberField(x^2 + 3)
            sage: K.trace_pairing([1,zeta3])
            [ 2  0]
            [ 0 -6]
        """
        import sage.matrix.matrix_space
        A = sage.matrix.matrix_space.MatrixSpace(self.base_ring(), len(v))(0)
        for i in range(len(v)):
            for j in range(i,len(v)):
                t = (self(v[i]*v[j])).trace()
                A[i,j] = t
                A[j,i] = t
        return A

    def uniformizer(self, P, others = "positive"):
        """
        Returns an element of self with valuation 1 at the prime ideal P.

        INPUT:
            self -- a number field
            P -- a prime ideal of self
            others -- either "positive" (default), in which case the element will have
                      non-negative valuation at all other primes of self,
                      or "negative", in which case the element will have non-positive
                      valuation at all other primes of self.

        NOTE: When P is principal (e.g. always when self has class number
        one) the result may or may not be a generator of P!

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 5); K
            Number Field in a with defining polynomial x^2 + 5
            sage: P,Q = K.ideal(3).prime_factors()
            sage: P
            Fractional ideal (3, a + 1)
            sage: pi=K.uniformizer(P); pi
            a + 1
            sage: K.ideal(pi).factor()
            (Fractional ideal (2, a + 1)) * (Fractional ideal (3, a + 1))
            sage: pi=K.uniformizer(P,'negative'); pi
            1/2*a + 1/2
            sage: K.ideal(pi).factor()
            (Fractional ideal (2, a + 1))^-1 * (Fractional ideal (3, a + 1))

            sage: K = CyclotomicField(9)
            sage: Plist=K.ideal(17).prime_factors()
            sage: pilist = [K.uniformizer(P) for P in Plist]
            sage: [pi.is_integral() for pi in pilist]
            [True, True, True]
            sage: [pi.valuation(P) for pi,P in zip(pilist,Plist)]
            [1, 1, 1]
            sage: [ pilist[i] in Plist[i] for i in range(len(Plist)) ]
            [True, True, True]
        """
        if not is_NumberFieldIdeal(P):
            P = self.ideal(P)
        if not P.is_maximal():
            raise ValueError, "P must be a nonzero prime"
        if others == "negative":
            P = ~P
        elif others != "positive":
            raise ValueError, "others must be 'positive' or 'negative'"
        nf = self.pari_nf()
        a = self(nf.idealappr(P.pari_hnf()))
        if others == "negative":
            a = ~a
        return a

    def units(self, proof=None):
        """
        Return generators for the unit group modulo torsion.

        ALGORITHM: Uses PARI's bnfunit command.

        INPUTS:
            proof -- default: True

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: A = x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3
            sage: K = NumberField(A, 'a')
            sage: K.units()
            [8/275*a^3 - 12/55*a^2 + 15/11*a - 2]

        Sage might not be able to provably compute the unit group:

            sage: K = NumberField(x^17 + 3, 'a')
            sage: K.units(proof=True) # default
            Traceback (most recent call last):
            ...
            PariError: not enough precomputed primes, need primelimit ~  (35)

        In this case, one can ask for the conjectural unit group (correct if
        the Generalized Riemann Hypothesis is true):

            sage: K.units(proof=False)
            [a^9 + a - 1,
            a^16 - a^15 + a^14 - a^12 + a^11 - a^10 - a^8 + a^7 - 2*a^6 + a^4 - 3*a^3 + 2*a^2 - 2*a + 1,
            2*a^16 - a^14 - a^13 + 3*a^12 - 2*a^10 + a^9 + 3*a^8 - 3*a^6 + 3*a^5 + 3*a^4 - 2*a^3 - 2*a^2 + 3*a + 4,
            2*a^16 - 3*a^15 + 3*a^14 - 3*a^13 + 3*a^12 - a^11 + a^9 - 3*a^8 + 4*a^7 - 5*a^6 + 6*a^5 - 4*a^4 + 3*a^3 - 2*a^2 - 2*a + 4,
            a^15 - a^12 + a^10 - a^9 - 2*a^8 + 3*a^7 + a^6 - 3*a^5 + a^4 + 4*a^3 - 3*a^2 - 2*a + 2,
            a^16 - a^15 - 3*a^14 - 4*a^13 - 4*a^12 - 3*a^11 - a^10 + 2*a^9 + 4*a^8 + 5*a^7 + 4*a^6 + 2*a^5 - 2*a^4 - 6*a^3 - 9*a^2 - 9*a - 7,
            a^15 + a^14 + 2*a^11 + a^10 - a^9 + a^8 + 2*a^7 - a^5 + 2*a^3 - a^2 - 3*a + 1,
            3*a^16 + 3*a^15 + 3*a^14 + 3*a^13 + 3*a^12 + 2*a^11 + 2*a^10 + 2*a^9 + a^8 - a^7 - 2*a^6 - 3*a^5 - 3*a^4 - 4*a^3 - 6*a^2 - 8*a - 8]

        The provable and the conjectural results are cached separately (this
        fixes trac \#2504):

            sage: K.units(proof=True)
            Traceback (most recent call last):
            ...
            PariError: not enough precomputed primes, need primelimit ~  (35)
        """
        proof = proof_flag(proof)

        # if we have cached provable results, return them immediately
        try:
            return self.__units
        except AttributeError:
            pass

        # if proof==False and we have cached results, return them immediately
        if not proof:
            try:
                return self.__units_no_proof
            except AttributeError:
                pass

        # get Pari to compute the units
        B = self.pari_bnf(proof).bnfunit()
        R = self.polynomial().parent()
        if proof:
            # cache the provable results and return them
            self.__units = [self(R(g)) for g in B]
            return self.__units
        else:
            # cache the conjectural results and return them
            self.__units_no_proof = [self(R(g)) for g in B]
            return self.__units_no_proof

    def zeta(self, n=2, all=False):
        """
        If all is False, return a primitive n-th root of unity in this
        field, or raise an ArithmeticError exception if there are none.

        If all is True, return a list of all primitive n-th roots of
        unity in this field (possibly empty).

        Note that if one wants to know the maximal root of unity
        in this field, one can use self.zeta_order().

        INPUT:
            n -- positive integer
            all -- bool, default: False.

        EXAMPLES:
            sage: K.<z> = NumberField(x^2 + 3)
            sage: K.zeta(1)
            1
            sage: K.zeta(2)
            -1
            sage: K.zeta(2, all=True)
            [-1]
            sage: K.zeta(3)
            1/2*z - 1/2
            sage: K.zeta(3, all=True)
            [1/2*z - 1/2, -1/2*z - 1/2]
            sage: K.zeta(4)
            Traceback (most recent call last):
            ...
            ArithmeticError: There are no 4-th roots of unity in self.

            sage: r.<x> = QQ[]
            sage: K.<b> = NumberField(x^2+1)
            sage: K.zeta(4)
            b
            sage: K.zeta(4,all=True)
            [b, -b]
            sage: K.zeta(3)
            Traceback (most recent call last):
            ...
            ArithmeticError: There are no 3-rd roots of unity in self.
            sage: K.zeta(3,all=True)
            []
        """
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n (=%s) must be positive"%n
        if n == 1:
            if all:
                return [self(1)]
            else:
                return self(1)
        elif n == 2:
            if all:
                return [self(-1)]
            else:
                return self(-1)
        else:
            field = self.absolute_field('a')
            from_field = field.structure()[0]
            f = field.polynomial_ring().cyclotomic_polynomial(n)
            F = field['x'](f)
            R = F.roots()
            if len(R) == 0:
                if all:
                    return []
                else:
                    if n == 1:
                        th = 'st'
                    elif n == 2:
                        th = 'nd'
                    elif n == 3:
                        th = 'rd'
                    else:
                        th = 'th'
                    raise ArithmeticError, "There are no %s-%s roots of unity in self."%(n,th)
            if all:
                return [from_field(r[0]) for r in R]
            else:
                return from_field(R[0][0])

    def zeta_order(self):
        r"""
        Return the number of roots of unity in this field.

        EXAMPLES:
            sage: F.<alpha> = NumberField(x**22+3)
            sage: F.zeta_order()
            6
            sage: F.<alpha> = NumberField(x**2-7)
            sage: F.zeta_order()
            2
        """
        return ZZ(self.pari_nf().nfrootsof1()[0])

    def number_of_roots_of_unity(self):
        """
        Return number of roots of unity in this field.

        EXAMPLES:
            sage: K.<b> = NumberField(x^2+1)
            sage: K.number_of_roots_of_unity()
            4
        """
        return ZZ(self.pari_nf().nfrootsof1()[0])

    def roots_of_unity(self):
        """
        Return all the roots of unity in this field, primitive or not.

        EXAMPLES:
            sage: K.<b> = NumberField(x^2+1)
            sage: zs = K.roots_of_unity(); zs
            [b, -1, -b, 1]
            sage: [ z**K.number_of_roots_of_unity() for z in zs ]
            [1, 1, 1, 1]
        """
        temp = self.pari_nf().nfrootsof1()
        n = ZZ(temp[0])
        z = self(temp[1])
        primitives = [ z**k for k in range(1, n) if sage.rings.arith.gcd(k,n)==1 ]
        primitives.sort(cmp=lambda z,w: len(str(z))-len(str(w)))
        z = primitives[0]
        return [ z**k for k in range(1, n+1) ]

    def zeta_coefficients(self, n):
        """
        Compute the first n coefficients of the Dedekind zeta function
        of this field as a Dirichlet series.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: NumberField(x^2+1, 'a').zeta_coefficients(10)
            [1, 1, 0, 1, 2, 0, 0, 1, 1, 2]
        """
        return self.pari_nf().dirzetak(n)


class NumberField_absolute(NumberField_generic):

    def __init__(self, polynomial, name, latex_name=None, check=True):
        """
        Function to initialize an absolute number field.
        """
        NumberField_generic.__init__(self, polynomial, name, latex_name, check)
        self._element_class = number_field_element.NumberFieldElement_absolute
        self._zero_element = self(0)
        self._one_element =  self(1)

    def _magma_convert_(self, magma):
        """
        Return Magma version of this number field.

        EXAMPLES:
            sage: R.<t> = PolynomialRing(RationalField())
            sage: K.<a> = NumberField(t^2 + 1)
            sage: L = magma(K)    # optional - magma
            sage: L               # optional - magma
            Number Field with defining polynomial t^2 + 1 over the Rational Field
            sage: L.1             # optional - magma
            a
            sage: L.1^2           # optional - magma
            -1
        """
        # Get magma version of defining polynomial of this number field
        f = magma(self.defining_polynomial())
        # Make number field from that polynomial
        K = f.NumberField()
        # Set variable name of the Magma version of this number field
        K.AssignNames(self.variable_names())
        return K

    def base_field(self):
        """
        Returns the base field of self, which is always QQ

        EXAMPLES:
            sage: K = CyclotomicField(5)
            sage: K.base_field()
            Rational Field
        """
        return QQ

    def is_absolute(self):
        """
        Returns True since self is an absolute field.

        EXAMPLES:
            sage: K = CyclotomicField(5)
            sage: K.is_absolute()
            True
        """
        return True

    def absolute_polynomial(self):
        r"""
        Return absolute polynomial that defines this absolute field.
        This is the same as \code{self.polynomial()}.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.absolute_polynomial ()
            x^2 + 1
        """
        return self.polynomial()

    def __reduce__(self):
        """
        TESTS:
            sage: Z = var('Z')
            sage: K.<w> = NumberField(Z^3 + Z + 1)
            sage: L = loads(dumps(K))
            sage: L
            Number Field in w with defining polynomial Z^3 + Z + 1
            sage: print L == K
            True
        """
        return NumberField_absolute_v1, (self.polynomial(), self.variable_name(), self.latex_variable_name())

    def optimized_representation(self, names=None, both_maps=True):
        """
        Return a field isomorphic to self with a better defining
        polynomial if possible, along with field isomorphisms from the
        new field to self and from self to the new field.

        EXAMPLES:
        We construct a compositum of 3 quadratic fields, then find an optimized
        representation and transform elements back and forth.
            sage: K = NumberField([x^2 + p for p in [5, 3, 2]],'a').absolute_field('b'); K
            Number Field in b with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
            sage: L, from_L, to_L = K.optimized_representation()
            sage: L    # your answer may different, since algorithm is random
            Number Field in a14 with defining polynomial x^8 + 4*x^6 + 7*x^4 + 36*x^2 + 81
            sage: to_L(K.0)   # random
            4/189*a14^7 - 1/63*a14^6 + 1/27*a14^5 + 2/9*a14^4 - 5/27*a14^3 + 8/9*a14^2 + 3/7*a14 + 3/7
            sage: from_L(L.0)   # random
            1/1152*a1^7 + 1/192*a1^6 + 23/576*a1^5 + 17/96*a1^4 + 37/72*a1^3 + 5/6*a1^2 + 55/24*a1 + 3/4

        The transformation maps are mutually inverse isomorphisms.
            sage: from_L(to_L(K.0))
            b
            sage: to_L(from_L(L.0))     # random
            a14
        """
        return self.optimized_subfields(degree=self.degree(), name=names, both_maps=both_maps)[0]

    def optimized_subfields(self, degree=0, name=None, both_maps=True):
        """
        Return optimized representations of many (but *not* necessarily
        all!)  subfields of self of degree 0, or of all possible
        degrees if degree is 0.

        EXAMPLES:
            sage: K = NumberField([x^2 + p for p in [5, 3, 2]],'a').absolute_field('b'); K
            Number Field in b with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
            sage: L = K.optimized_subfields(name='b')
            sage: L[0][0]
            Number Field in b0 with defining polynomial x - 1
            sage: L[1][0]
            Number Field in b1 with defining polynomial x^2 - x + 1
            sage: [z[0] for z in L]          # random -- since algorithm is random
            [Number Field in b0 with defining polynomial x - 1,
             Number Field in b1 with defining polynomial x^2 - x + 1,
             Number Field in b2 with defining polynomial x^4 - 5*x^2 + 25,
             Number Field in b3 with defining polynomial x^4 - 2*x^2 + 4,
             Number Field in b4 with defining polynomial x^8 + 4*x^6 + 7*x^4 + 36*x^2 + 81]

        We examine one of the optimized subfields in more detail:
             sage: M, from_M, to_M = L[2]
             sage: M                             # random
             Number Field in b2 with defining polynomial x^4 - 5*x^2 + 25
             sage: from_M     # may be slightly random
             Ring morphism:
               From: Number Field in b2 with defining polynomial x^4 - 5*x^2 + 25
               To:   Number Field in a1 with defining polynomial x^8 + 40*x^6 + 352*x^4 + 960*x^2 + 576
               Defn: b2 |--> -5/1152*a1^7 + 1/96*a1^6 - 97/576*a1^5 + 17/48*a1^4 - 95/72*a1^3 + 17/12*a1^2 - 53/24*a1 - 1

        The to_M map is None, since there is no map from K to M:
             sage: to_M

        We apply the from_M map to the generator of M, which gives a rather
        large element of $K$:
             sage: from_M(M.0)          # random
             -5/1152*a1^7 + 1/96*a1^6 - 97/576*a1^5 + 17/48*a1^4 - 95/72*a1^3 + 17/12*a1^2 - 53/24*a1 - 1

        Nevertheless, that large-ish element lies in a degree 4 subfield:
             sage: from_M(M.0).minpoly()   # random
             x^4 - 5*x^2 + 25
        """
        return self._subfields_helper(degree=degree,name=name,
                                      both_maps=both_maps,optimize=True)

    def change_names(self, names):
        r"""
        Return number field isomorphic to self but with the given
        generator name.

        INPUT:
            names -- should be exactly one variable name.

        Also, \code{K.structure()} returns from_K and to_K, where
        from_K is an isomorphism from K to self and to_K is an
        isomorphism from self to K.

        EXAMPLES:
            sage: K.<z> = NumberField(x^2 + 3); K
            Number Field in z with defining polynomial x^2 + 3
            sage: L.<ww> = K.change_names()
            sage: L
            Number Field in ww with defining polynomial x^2 + 3
            sage: L.structure()[0]
            Number field isomorphism from Number Field in ww with defining polynomial x^2 + 3 to Number Field in z with defining polynomial x^2 + 3 given by variable name change
            sage: L.structure()[0](ww + 5/3)
            z + 5/3
        """
        return self.absolute_field(names)

    def subfields(self, degree=0, name=None):
        """
        EXAMPLES:
            sage: K.<a> = NumberField( [x^3 - 2, x^2 + x + 1] )
            sage: K = K.absolute_field('b')
            sage: S = K.subfields()
            sage: len(S)
            6
            sage: [k[0].polynomial() for k in S]
            [x - 3,
             x^2 - 3*x + 9,
             x^3 - 3*x^2 + 3*x + 1,
             x^3 - 3*x^2 + 3*x + 1,
             x^3 - 3*x^2 + 3*x - 17,
             x^6 - 3*x^5 + 6*x^4 - 11*x^3 + 12*x^2 + 3*x + 1]
             sage: R.<t> = QQ[]
             sage: L = NumberField(t^3 - 3*t + 1, 'c')
             sage: [k[1] for k in L.subfields()]
             [Ring morphism:
               From: Number Field in c0 with defining polynomial t
               To:   Number Field in c with defining polynomial t^3 - 3*t + 1
               Defn: 0 |--> 0,
              Ring morphism:
               From: Number Field in c1 with defining polynomial t^3 - 3*t + 1
               To:   Number Field in c with defining polynomial t^3 - 3*t + 1
               Defn: c1 |--> c]
        """
        return self._subfields_helper(degree=degree, name=name,
                                      both_maps=True, optimize=False)

    def _subfields_helper(self, degree=0, name=None, both_maps=True, optimize=False):
        """
        Internal function: common code for optimized_subfields() and subfields().
        """
        if name is None:
            name = self.variable_names()
        name = sage.structure.parent_gens.normalize_names(1, name)[0]
        try:
            return self.__subfields[name, degree, both_maps, optimize]
        except AttributeError:
            self.__subfields = {}
        except KeyError:
            pass
        f = pari(self.polynomial())
        if optimize:
            v = f.polred(2)
            elts = v[0]; polys = v[1]
        else:
            v = f.nfsubfields(degree)
            elts = [x[1] for x in v]; polys = [x[0] for x in v]

        R = self.polynomial_ring()

        ans = []
        for i in range(len(elts)):
            f = R(polys[i])
            if not (degree == 0 or f.degree() == degree):
                continue
            a = self(R(elts[i]))
            K = NumberField(f, names=name + str(i))

            from_K = K.hom([a])    # check=False here ??   would be safe unless there are bugs.

            if both_maps and K.degree() == self.degree():
                g = K['x'](self.polynomial())
                v = g.roots()
                a = from_K(K.gen())
                for i in range(len(v)):
                    r = g.roots()[i][0]
                    to_K = self.hom([r])    # check=False here ??
                    if to_K(a) == K.gen():
                        break
            else:
                to_K = None
            ans.append((K, from_K, to_K))
        ans = Sequence(ans, immutable=True, cr=True)
        self.__subfields[name, degree, both_maps, optimize] = ans
        return ans


    def maximal_order(self, v=None):
        """
        Return the maximal order, i.e., the ring of integers, associated
        to this number field.

        INPUT:
            v -- (default: None) None, a prime, or a list of primes.
                 * if v is None, return the maximal order.
                 * if v is a prime, return an order that is p-maximal.
                 * if v is a list, return an order that is maximal at
                   each prime in the list v.

        EXAMPLES:
        In this example, the maximal order cannot be generated
        by a single element.
            sage: k.<a> = NumberField(x^3 + x^2 - 2*x+8)
            sage: o = k.maximal_order()
            sage: o
            Maximal Order in Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8

        We compute $p$-maximal orders for several $p$.  Note that computing
        a $p$-maximal order is much faster in general than computing
        the maximal order:
            sage: p = next_prime(10^22); q = next_prime(10^23)
            sage: K.<a> = NumberField(x^3 - p*q)
            sage: K.maximal_order([3]).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]
            sage: K.maximal_order([2]).basis()
            [1, a, a^2]
            sage: K.maximal_order([p]).basis()
            [1, a, a^2]
            sage: K.maximal_order([q]).basis()
            [1, a, a^2]
            sage: K.maximal_order([p,3]).basis()
            [1/3*a^2 + 1/3*a + 1/3, a, a^2]

        An example with bigger discriminant:
            sage: p = next_prime(10^97); q = next_prime(10^99)
            sage: K.<a> = NumberField(x^3 - p*q)
            sage: K.maximal_order(prime_range(10000)).basis()
            [1, a, a^2]
        """
        v = self._normalize_prime_list(v)

        try:
            return self.__maximal_order[v]
        except AttributeError:
            self.__maximal_order = {}
        except KeyError:
            pass

        B = self._compute_integral_basis(v = v)

        if len(v) == 0 or v is None:
            is_maximal = True
        else:
            is_maximal = False

        import sage.rings.number_field.order as order
        O = order.absolute_order_from_module_generators(B,
                 check_integral=False, check_rank=False,
                 check_is_ring=False, is_maximal = is_maximal)

        self.__maximal_order[v] = O
        return O

    def order(self, *gens, **kwds):
        r"""
        Return the order with given ring generators in the maximal
        order of this number field.

        INPUT:
            gens -- list of elements of self; if no generators are
                    given, just returns the cardinality of this number
                    field (oo) for consistency.
            check_is_integral -- bool (default: True), whether to check
                  that each generator is integral.
            check_rank -- bool (default: True), whether to check that
                  the ring generated by gens is of full rank.
            allow_subfield -- bool (default: False), if True and the generators
                  do not generate an order, i.e., they generate a subring
                  of smaller rank, instead of raising an error, return
                  an order in a smaller number field.

        EXAMPLES:
            sage: k.<i> = NumberField(x^2 + 1)
            sage: k.order(2*i)
            Order in Number Field in i with defining polynomial x^2 + 1
            sage: k.order(10*i)
            Order in Number Field in i with defining polynomial x^2 + 1
            sage: k.order(3)
            Traceback (most recent call last):
            ...
            ValueError: the rank of the span of gens is wrong
            sage: k.order(i/2)
            Traceback (most recent call last):
            ...
            ValueError: each generator must be integral

        Alternatively, an order can be constructed by adjoining
        elements to $\ZZ$:

        """
        if len(gens) == 0:
            return NumberField_generic.order(self)
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        gens = [self(x) for x in gens]
        import sage.rings.number_field.order as order
        return order.absolute_order_from_ring_generators(gens, **kwds)

    def vector_space(self):
        """
        Return a vector space V and isomorphisms self --> V and V --> self.

        OUTPUT:
            V -- a vector space over the rational numbers
            from_V -- an isomorphism from V to self
            to_V -- an isomorphism from self to V

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 + 2)
            sage: V, from_V, to_V  = k.vector_space()
            sage: from_V(V([1,2,3]))
            3*a^2 + 2*a + 1
            sage: to_V(1 + 2*a + 3*a^2)
            (1, 2, 3)
            sage: V
            Vector space of dimension 3 over Rational Field
            sage: to_V
            Isomorphism from Number Field in a with defining polynomial x^3 + 2 to Vector space of dimension 3 over Rational Field
            sage: from_V(to_V(2/3*a - 5/8))
            2/3*a - 5/8
            sage: to_V(from_V(V([0,-1/7,0])))
            (0, -1/7, 0)
        """
        try:
            return self.__vector_space
        except AttributeError:
            V = QQ**self.degree()
            from_V = maps.MapVectorSpaceToNumberField(V, self)
            to_V   = maps.MapNumberFieldToVectorSpace(self, V)
            self.__vector_space = (V, from_V, to_V)
            return self.__vector_space

    def absolute_vector_space(self):
        r"""
        Return vector space over $\QQ$ corresponding to this number
        field, along with maps from that space to this number field
        and in the other direction.

        For an absolute extension this is identical to
        \code{self.vector_space()}.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 - 5)
            sage: K.absolute_vector_space()
            (Vector space of dimension 3 over Rational Field,
             Isomorphism from Vector space of dimension 3 over Rational Field to Number Field in a with defining polynomial x^3 - 5,
             Isomorphism from Number Field in a with defining polynomial x^3 - 5 to Vector space of dimension 3 over Rational Field)
        """
        return self.vector_space()

    def galois_closure(self, names=None):
        """
        Return number field $K$ that is the Galois closure of self,
        i.e., is generated by all roots of the defining polynomial of
        self

        INPUT:
            names -- variable name for Galois closure

        EXAMPLES:
            sage: K.<a> = NumberField(x^4 - 2)
            sage: M = K.galois_closure('b'); M
            Number Field in b with defining polynomial x^8 + 28*x^4 + 2500
            sage: L.<a2> = K.galois_closure(); L
            Number Field in a2 with defining polynomial x^8 + 28*x^4 + 2500
            sage: K.galois_group().order()
            8

            sage: phi = K.embeddings(L)[0]
            sage: phi(K.0)
            1/120*a2^5 + 19/60*a2
            sage: phi(K.0).minpoly()
            x^4 - 2
        """
        try:
            return self.__galois_closure.change_names(names)
        except AttributeError:
            pass
        G = self.galois_group()
        K = self
        while K.degree() < G.order():
            K = K.composite_fields(self, names=names)[-1]
        self.__galois_closure = K.change_names(names)
        return self.__galois_closure

    def embeddings(self, K):
        """
        Compute all field embeddings of self into the field K (which
        need not even be a number field, e.g., it could be the complex
        numbers). This will return an identical result when given K as
        input again.

        If possible, the most natural embedding of K into self
        is put first in the list.

        INPUT:
            K -- a number field

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 - 2)
            sage: L.<a1> = K.galois_closure(); L
            Number Field in a1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: K.embeddings(L)[0]
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Number Field in a1 with defining polynomial x^6 + 40*x^3 + 1372
              Defn: a |--> 1/84*a1^4 + 13/42*a1
            sage: K.embeddings(L)  is K.embeddings(L)
            True

        We embed a quadratic field into a cyclotomic field:
            sage: L.<a> = QuadraticField(-7)
            sage: K = CyclotomicField(7)
            sage: L.embeddings(K)
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^2 + 7
              To:   Cyclotomic Field of order 7 and degree 6
              Defn: a |--> 2*zeta7^4 + 2*zeta7^2 + 2*zeta7 + 1,
            Ring morphism:
              From: Number Field in a with defining polynomial x^2 + 7
              To:   Cyclotomic Field of order 7 and degree 6
              Defn: a |--> -2*zeta7^4 - 2*zeta7^2 - 2*zeta7 - 1
            ]

        We embed a cubic field in the complex numbers:
            sage: K.<a> = NumberField(x^3 - 2)
            sage: K.embeddings(CC)
            [
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Complex Field with 53 bits of precision
              Defn: a |--> -0.62996052494743... - 1.09112363597172*I,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Complex Field with 53 bits of precision
              Defn: a |--> -0.62996052494743... + 1.09112363597172*I,
            Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Complex Field with 53 bits of precision
              Defn: a |--> 1.25992104989487
            ]
        """
        try:
            return self.__embeddings[K]
        except AttributeError:
            self.__embeddings = {}
        except KeyError:
            pass
        f = K['x'](self.defining_polynomial())
        r = f.roots(); r.sort()
        v = [self.hom([e[0]], check=False) for e in r]
        # If there is an embedding that preserves variable names
        # then it is most natural, so we put it first.
        put_natural_embedding_first(v)

        self.__embeddings[K] = Sequence(v, cr=v!=[], immutable=True,
                                        check=False, universe=self.Hom(K))
        return self.__embeddings[K]

    def Minkowski_embedding(self, B=None, prec=None):
        r"""
        Return an nxn matrix over RDF whose columns are the images of
        the basis $\{1, \alpha, \dots, \alpha^{n-1}\}$ of self over
        $\QQ$ (as vector spaces), where here $\alpha$ is the
        generator of self over $\QQ$, i.e.  self.gen(0).  If B
        is not None, return the images of the vectors in B as the
        columns instead. If prec is not None, use RealField(prec)
        instead of RDF.

        This embedding is the so-called "Minkowski embedding" of a
        number field in $\RR^n$: given the $n$ embeddings
        $\sigma_1, \dots, \sigma_n$ of self in $\CC$, write
        $\sigma_1, \dots, \sigma_r$ for the real embeddings, and
        $\sigma_{r+1}, \dots, \sigma_{r+s}$ for choices of one of each
        pair of complex conjugate embeddings (in our case, we simply
        choose the one where the image of $\alpha$ has positive real
        part). Here $(r,s)$ is the signature of self.  Then the
        Minkowski embedding is given by:

          x |--> ( $\sigma_1(x)$, $\dots$, $\sigma_r(x)$,
                   $\sqrt{2}\Re(\sigma_{r+1}(x))$,
                   $\sqrt{2}\Im(\sigma_{r+1}(x))$,
                   $\dots$,
                   $\sqrt{2}\Re(\sigma_{r+s}(x))$,
                   $\sqrt{2}\Im(\sigma_{r+s}(x))$)

        Equivalently, this is an embedding of self in $\RR^n$
        so that the usual norm on $\RR^n$ coincides with
          $\|x\| = \sum_i |\sigma_i(x)|^2$
        on self.

        TODO: This could be much improved by implementing
        homomorphisms over VectorSpaces.

        EXAMPLES:
            sage: F.<alpha> = NumberField(x^3+2)
            sage: F.Minkowski_embedding()
            [ 1.00000000000000 -1.25992104989487  1.58740105196820]
            [ 1.41421356237... 0.8908987181... -1.12246204830...]
            [0.000000000000000  1.54308184421...  1.94416129723...]
            sage: F.Minkowski_embedding([1, alpha+2, alpha^2-alpha])
            [ 1.00000000000000 0.740078950105127  2.84732210186307]
            [ 1.41421356237...  3.7193258428... -2.01336076644...]
            [0.000000000000000  1.54308184421... 0.40107945302...]
            sage: F.Minkowski_embedding() * (alpha + 2).vector().transpose()
            [0.740078950105127]
            [ 3.7193258428...]
            [ 1.54308184421...]
        """
        n = self.degree()
        if prec is None:
            R = sage.rings.real_double.RDF
        else:
            R = sage.rings.real_mpfr.RealField(prec)
        r,s = self.signature()
        places = self.places(prec=prec)

        if B is None:
            B = [ (self.gen(0))**i for i in range(n) ]

        A = ZZ['x']
        f = A.gen(0)**2-2
        sqrt2 = f.roots(R)[1][0]

        d = {}

        for col in range(n):

            for row in range(r):
                d[(row,col)] = places[row](B[col])

            for i in range(s):
                z = places[r+i](B[col])
                d[(r+2*i,col)] = z.real()*sqrt2
                d[(r+2*i+1,col)] = z.imag()*sqrt2


        M = sage.matrix.all.matrix(d)

        return M


    def places(self, all_complex=False, prec=None):
        """
        Return the collection of all places of self. By default, this
        returns the set of real places as homomorphisms into RIF
        first, followed by a choice of one of each pair of complex
        conjugate homomorphisms into CIF.

        On the other hand, if prec is not None, we simply return
        places into RealField(prec) and ComplexField(prec) (or RDF,
        CDF if prec=53).

        There is an optional flag all_complex, which defaults to
        False. If all_complex is True, then the real embeddings are
        returned as embeddings into CIF instead of RIF.

        EXAMPLES:
            sage: F.<alpha> = NumberField(x^3-100*x+1) ; F.places()
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 - 100*x + 1
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> -10.00499625499181184573367219280,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 - 100*x + 1
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> 0.01000001000003000012000055000273,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 - 100*x + 1
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> 9.994996244991781845613530439509]

            sage: F.<alpha> = NumberField(x^3+7) ; F.places()
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> -1.912931182772389101199116839549,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 53 bits of precision
            Defn: alpha |--> 0.956465591386195 + 1.65664699997230*I]

            sage: F.<alpha> = NumberField(x^3+7) ; F.places(all_complex=True)
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 53 bits of precision
            Defn: alpha |--> -1.91293118277239,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 53 bits of precision
            Defn: alpha |--> 0.956465591386195 + 1.65664699997230*I]
            sage: F.places(prec=10)
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Real Field with 10 bits of precision
            Defn: alpha |--> -1.9,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^3 + 7
            To:   Complex Field with 10 bits of precision
            Defn: alpha |--> 0.96 + 1.7*I]
        """
        if prec is None:
            R = RIF
            C = CIF
        elif prec == 53:
            R = sage.rings.real_double.RDF
            C = sage.rings.complex_double.CDF
        else:
            R = sage.rings.real_mpfr.RealField(prec)
            C = sage.rings.complex_field.ComplexField(prec)

        ## first, find the intervals with roots, and see how much
        ## precision we need to approximate the roots
        ##
        all_intervals = [ x[0] for x in self.defining_polynomial().roots(C) ]

        ## first, set up the real places
        if all_complex:
            real_intervals = [ x for x in all_intervals if x.imag().is_zero() ]
        else:
            real_intervals = [ x[0] for x in self.defining_polynomial().roots(R) ]

        if prec is None:
            real_places = [ self.hom([i.center()], check=False) for i in real_intervals ]

            complex_places = [ self.hom([i.center()], check=False) for i in
                               all_intervals if i.imag() > 0 ]
        else:
            real_places = [ self.hom([i], check=False) for i in real_intervals ]

            complex_places = [ self.hom([i], check=False) for i in
                               all_intervals if i.imag() > 0 ]

        return real_places + complex_places

    def real_places(self, prec=None):
        """
        Return all real places of self as homomorphisms into RIF.

        EXAMPLES:
            sage: F.<alpha> = NumberField(x^4-7) ; F.real_places()
            [Ring morphism:
            From: Number Field in alpha with defining polynomial x^4 - 7
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> -1.626576561697785743211232345494,
            Ring morphism:
            From: Number Field in alpha with defining polynomial x^4 - 7
            To:   Real Field with 106 bits of precision
            Defn: alpha |--> 1.626576561697785743211232345494]
        """
        return self.places(prec=prec)[0:self.signature()[0]]


    def relativize(self, alpha, names):
        r"""
        Given an element alpha in self, return a relative number field
        $K$ isomorphic to self that is relative over the absolute field
        $\QQ(\alpha)$, along with isomorphisms from $K$ to self and
        from self to K.

        INPUT:
            alpha -- an element of self.
            names -- 2-tuple of names of generator for output
                     field K and the subfield QQ(alpha)
                     names[0] generators K and names[1] QQ(alpha).

        OUTPUT:
            K   -- relative number field

        Also, \code{K.structure()} returns from_K and to_K, where
        from_K is an isomorphism from K to self and to_K is an isomorphism
        from self to K.

        EXAMPLES:
            sage: K.<a> = NumberField(x^10 - 2)
            sage: L.<c,d> = K.relativize(a^4 + a^2 + 2); L
            Number Field in c with defining polynomial x^2 - 1/5*d^4 + 8/5*d^3 - 23/5*d^2 + 7*d - 18/5 over its base field
            sage: c.absolute_minpoly()
            x^10 - 2
            sage: d.absolute_minpoly()
            x^5 - 10*x^4 + 40*x^3 - 90*x^2 + 110*x - 58
            sage: (a^4 + a^2 + 2).minpoly()
            x^5 - 10*x^4 + 40*x^3 - 90*x^2 + 110*x - 58
            sage: from_L, to_L = L.structure()
            sage: to_L(a)
            c
            sage: to_L(a^4 + a^2 + 2)
            d
            sage: from_L(to_L(a^4 + a^2 + 2))
            a^4 + a^2 + 2
        """
        # step 1: construct the abstract field generated by alpha.
        # step 2: make a relative extension of it.
        # step 3: construct isomorphisms

        names = sage.structure.parent_gens.normalize_names(2, names)

        # make sure alpha is in self
        alpha = self(alpha)

        f = alpha.minpoly()
        L = NumberField(f, names[1])

        g = self.defining_polynomial()
        h = L['x'](g)
        F = h.factor()

        for f, e in F:
            if L.degree() * f.degree() == self.degree():
                M = L.extension(f, names[0])
                beta = M(L.gen())
                try:
                    to_M = self.hom([M.gen(0)], M, check=True)  # be paranoid
                except TypeError:
                    continue
                if to_M(alpha) == beta:
                    # Bingo.
                    # We have now constructed a relative
                    # number field M, and an isomorphism
                    # self --> M that sends alpha to
                    # the generator of the intermediate field.
                    from_M = M.hom([self.gen()], self, check=True)
                    M._set_structure(from_M, to_M)  # don't have to
                                                    # worry about caching since relative number fields aren't cached.
                    return M

        assert False, "bug in relativize"









class NumberField_relative(NumberField_generic):
    """
    EXAMPLES:
        sage: K.<a> = NumberField(x^3 - 2)
        sage: t = K['x'].gen()
        sage: L.<b> = K.extension(t^2+t+a); L
        Number Field in b with defining polynomial x^2 + x + a over its base field
    """
    def __init__(self, base, polynomial, name,
                 latex_name=None, names=None, check=True):
        r"""
        INPUT:
            base -- the base field
            polynomial -- must be defined in the ring \code{K['x']}, where
                          K is the base field.
            name -- variable name
            latex_name -- latex variable name
            names --
            check -- whether to check irreducibility of polynomial.

        EXAMPLES:
            sage: K.<x> = CyclotomicField(5)[]
            sage: W.<a> = NumberField(x^2 + 1)
            sage: W
            Number Field in a with defining polynomial x^2 + 1 over its base field
            sage: type(W)
            <class 'sage.rings.number_field.number_field.NumberField_relative'>

        Test that check=False really skips the test:
            sage: W.<a> = NumberField(K.cyclotomic_polynomial(5), check=False)
            sage: W
            Number Field in a with defining polynomial x^4 + x^3 + x^2 + x + 1 over its base field

        A relative extension of a relative extension:
            sage: x = var('x')
            sage: k.<a> = NumberField([x^2 + 2, x^2 + 1])
            sage: l.<b> = k.extension(x^2 + 3)
            sage: l
            Number Field in b with defining polynomial x^2 + 3 over its base field
            sage: l.base_field()
            Number Field in a0 with defining polynomial x^2 + 2 over its base field
            sage: l.base_field().base_field()
            Number Field in a1 with defining polynomial x^2 + 1
        """
        if not names is None: name = names
        if not is_NumberField(base):
            raise TypeError, "base (=%s) must be a number field"%base
        if not isinstance(polynomial, polynomial_element.Polynomial):
            try:
                polynomial = polynomial.polynomial(base)
            except (AttributeError, TypeError), msg:
                raise TypeError, "polynomial (=%s) must be a polynomial."%repr(polynomial)
        if name == base.variable_name():
            raise ValueError, "Base field and extension cannot have the same name"
        if polynomial.parent().base_ring() != base:
            polynomial = polynomial.change_ring(base)
            #raise ValueError, "The polynomial must be defined over the base field"

        # Generate the nf and bnf corresponding to the base field
        # defined as polynomials in y, e.g. for rnfisfree

        # Convert the polynomial defining the base field into a
        # polynomial in y to satisfy PARI's ordering requirements.

        if base.is_relative():
            abs_base = base.absolute_field('a')
            from_abs_base, to_abs_base = abs_base.structure()
        else:
            abs_base = base
            from_abs_base = maps.IdentityMap(base)
            to_abs_base = maps.IdentityMap(base)

        self.__absolute_base_field = abs_base, from_abs_base, to_abs_base
        Qx = abs_base.polynomial().parent()
        Qy = (abs_base.polynomial().base_ring())['y']
        phi = Qx.hom([Qy.gen()])
        base_polynomial_y = phi(abs_base.polynomial())

        self.__base_nf = pari(base_polynomial_y).nfinit()
        self.__base_bnf = pari(base_polynomial_y).bnfinit()

        # Use similar methods to convert the polynomial defining the
        # relative extension into a polynomial in x, with y denoting
        # the generator of the base field.
        # NOTE: This should be rewritten if there is a way to extend
        #       homomorphisms K -> K' to homomorphisms K[x] -> K'[x].

        base_field_y = NumberField(abs_base.polynomial(), 'y')
        Kx = base_field_y['x']
        i = abs_base.hom([base_field_y.gen()]) # inclusion K -> K' with a -> y
        rel_coeffs = [i(to_abs_base(c)) for c in polynomial.coeffs()]
        polynomial_y = Kx(rel_coeffs)

        if check:
            if not polynomial_y.is_irreducible():
                raise ValueError, "defining polynomial (%s) must be irreducible"%polynomial


        self.__pari_relative_polynomial = pari(str(polynomial_y))
        self.__rnf = self.__base_nf.rnfinit(self.__pari_relative_polynomial)

        self.__base_field = base
        self.__relative_polynomial = polynomial
        self.__pari_bnf_certified = False
        self._element_class = number_field_element.NumberFieldElement_relative

        self.__gens = [None]

        v = [None]
        K = base
        names = [name]
        while K != QQ:
            names.append(K.variable_name())
            v.append(K.gen())
            K = K.base_field()

        self._assign_names(tuple(names), normalize=False)

        NumberField_generic.__init__(self, self.absolute_polynomial(), name=None,
                                     latex_name=latex_name, check=False)

        v[0] = self._gen_relative()
        v = [self(x) for x in v]
        self.__gens = tuple(v)
        self._zero_element = self(0)
        self._one_element =  self(1)

    def change_names(self, names):
        r"""
        Return relative number field isomorphic to self but with the
        given generator names.

        INPUT:
            names -- number of names should be at most the number of
                     generators of self, i.e., the number of steps in
                     the tower of relative fields.

        Also, \code{K.structure()} returns from_K and to_K, where
        from_K is an isomorphism from K to self and to_K is an
        isomorphism from self to K.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: L.<c,d> = K.change_names()
            sage: L
            Number Field in c with defining polynomial x^4 + 3 over its base field
            sage: L.base_field()
            Number Field in d with defining polynomial x^2 + 2

        An example with a 3-level tower:
            sage: K.<a,b,c> = NumberField([x^2 + 17, x^2 + x + 1, x^3 - 2]); K
            Number Field in a with defining polynomial x^2 + 17 over its base field
            sage: L.<m,n,r> = K.change_names()
            sage: L
            Number Field in m with defining polynomial x^2 + 17 over its base field
            sage: L.base_field()
            Number Field in n with defining polynomial x^2 + x + 1 over its base field
            sage: L.base_field().base_field()
            Number Field in r with defining polynomial x^3 - 2
        """
        if len(names) == 0:
            names = self.variable_names()
        elif isinstance(names, str):
            names = names.split(',')
        K = self.base_field().change_names(tuple(names[1:]))
        L = K.extension(self.defining_polynomial(), names=names[0])
        return L

    def is_absolute(self):
        """
        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: K.is_absolute()
            False
            sage: K.is_relative()
            True
        """
        return False

    def gens(self):
        """
        Return the generators of this relative number field.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: K.gens()
            (a, b)
        """
        return self.__gens

    def ngens(self):
        """
        Return the number of generators of this relative number field.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: K.gens()
            (a, b)
            sage: K.ngens()
            2
        """
        return len(self.__gens)

    def gen(self, n=0):
        """
        Return the n'th generator of this relative number field.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: K.gens()
            (a, b)
            sage: K.gen(0)
            a
        """
        if n < 0 or n >= len(self.__gens):
            raise IndexError, "invalid generator %s"%n
        return self.__gens[n]

    def galois_closure(self, names=None):
        """
        Return the absolute number field $K$ that is the Galois
        closure of this relative number field.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: K.galois_closure('c')
            Number Field in c with defining polynomial x^16 + 144*x^14 + 8988*x^12 + 329616*x^10 + 7824006*x^8 + 113989680*x^6 + 1360354716*x^4 + 3470308272*x^2 + 9407642049
        """
        return self.absolute_field('a').galois_closure(names=names)

    def absolute_degree(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField([x^2 + 3, x^2 + 2])
            sage: K.absolute_degree()
            4
            sage: K.degree()
            2
        """
        return self.absolute_polynomial().degree()

    def maximal_order(self):
        """
        Return the maximal order, i.e., the ring of integers of this
        number field.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order(); OK.basis()
            [1, 1/2*a - 1/2*b, -1/2*b*a + 1/2, a]
            sage: charpoly(OK.1)
            x^2 + b*x + 1
            sage: charpoly(OK.2)
            x^2 - x + 1
            sage: O2 = K.order([3*a, 2*b])
            sage: O2.index_in(OK)
            144
        """
        try:
            return self.__maximal_order
        except AttributeError:
            pass
        K = self.absolute_field('a')
        from_K,_ = K.structure()
        O = K.maximal_order()
        B = [from_K(z) for z in O.basis()]
        OK = self.order(B, check_is_integral=False, check_rank=False)
        self.__maximal_order = OK
        return OK


    def __reduce__(self):
        """
        TESTS:
            sage: Z = var('Z')
            sage: K.<w> = NumberField(Z^3 + Z + 1)
            sage: L.<z> = K.extension(Z^3 + 2)
            sage: L = loads(dumps(K))
            sage: print L
            Number Field in w with defining polynomial Z^3 + Z + 1
            sage: print L == K
            True
        """
        return NumberField_relative_v1, (self.__base_field, self.polynomial(), self.variable_name(),
                                          self.latex_variable_name())

    def _repr_(self):
        """
        Return string representation of this relative number field.

        The base field is not part of the string representation.  To
        find out what the base field is use \code{self.base_field()}.

        EXAMPLES:
            sage: k.<a, b> = NumberField([x^5 + 2, x^7 + 3])
            sage: k
            Number Field in a with defining polynomial x^5 + 2 over its base field
            sage: k.base_field()
            Number Field in b with defining polynomial x^7 + 3
        """

        return "Number Field in %s with defining polynomial %s over its base field"%(self.variable_name(), self.polynomial())

        #return "Extension by %s of the Number Field in %s with defining polynomial %s"%(
        #self.polynomial(), self.base_field().variable_name(),
        #    self.base_field().polynomial())

    def _Hom_(self, codomain, cat=None):
        """
        Return homset of homomorphisms from this relative number field
        to the codomain.

        The cat option is currently ignored.   The result is not cached.

        EXAMPLES:
        This function is implicitly called by the Hom method or function.
            sage: K.<a,b> = NumberField([x^3 - 2, x^2+1])
            sage: K.Hom(K)
            Automorphism group of Number Field in a with defining polynomial x^3 - 2 over its base field
            sage: type(K.Hom(K))
            <class 'sage.rings.number_field.morphism.RelativeNumberFieldHomset'>
        """
        import morphism
        return morphism.RelativeNumberFieldHomset(self, codomain)

    def _latex_(self):
        r"""
        Return a \LaTeX representation of the extension.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^3 - 2)
            sage: t = K['x'].gen()
            sage: K.extension(t^2+t+a, 'b')._latex_()
            '( \\mathbf{Q}[a]/(a^{3} - 2) )[b]/(b^{2} + b + a)'
        """
        return "( %s )[%s]/(%s)"%(latex(self.base_field()), self.latex_variable_name(),
                              self.polynomial()._latex_(self.latex_variable_name()))

    def __call__(self, x):
        """
        Coerce x into this relative number field.

        EXAMPLES:
        We construct the composite of three quadratic fields, then
        coerce from the quartic subfield of the relative extension:

            sage: k.<a,b,c> = NumberField([x^2 + 5, x^2 + 3, x^2 + 1])
            sage: m = k.base_field(); m
            Number Field in b with defining polynomial x^2 + 3 over its base field
            sage: k(m.0)
            b
            sage: k(2/3)
            2/3
            sage: k(m.0^4)
            9

        TESTS:
            sage: K.<a> = NumberField(ZZ['x'].0^2 + 2, 'a')
            sage: L.<b> = K.extension(ZZ['x'].0 - a, 'b')
            sage: L(a)
            a
            sage: L(b+a)
            2*a
            sage: K.<a> = NumberField(ZZ['x'].0^5 + 2, 'a')
            sage: L.<b> = K.extension(ZZ['x'].0 - a, 'b')
            sage: L(a)
            a
            sage: L(a**3)
            a^3
            sage: L(a**2+b)
            a^2 + a
            sage: L.<b> = K.extension(ZZ['x'].0 + a/2, 'b')
            sage: L(a)
            a
            sage: L(b)
            -1/2*a
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            P = x.parent()
            from sage.rings.number_field.order import is_NumberFieldOrder
            if P is self:
                return x
            elif is_NumberFieldOrder(P) and P.number_field() is self:
                return self._element_class(self, x.polynomial())
            elif P == self:
                return self._element_class(self, x.polynomial())
            return self.__base_inclusion(self.base_field()(x))

        if not isinstance(x, (int, long, rational.Rational,
                              integer.Integer, pari_gen,
                              polynomial_element.Polynomial,
                              list)):
            return self.base_field()(x)

        return self._element_class(self, x)

    def _coerce_impl(self, x):
        """
        Canonical implicit coercion of x into self.

        Elements of this field canonically coerce in, as does anything
        that coerces into the base field of this field.

        EXAMPLES:
            sage: k.<a> = NumberField([x^5 + 2, x^7 + 3])
            sage: b = k(k.base_field().gen())
            sage: b = k._coerce_impl(k.base_field().gen())
            sage: b^7
            -3
            sage: k._coerce_impl(2/3)
            2/3
            sage: c = a + b  # this works
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            from sage.rings.number_field.order import is_NumberFieldOrder
            if x.parent() is self:
                return x
            elif is_NumberFieldOrder(x.parent()) and x.parent().number_field() is self:
                return self._element_class(self, x.polynomial())
            else:
                return self.__base_inclusion(x)

        return self.__base_inclusion(self.base_field()._coerce_impl(x))

    def __base_inclusion(self, element):
        """
        Given an element of the base field, give its inclusion into
        this extension in terms of the generator of this field.

        This is called by the canonical coercion map on elements from
        the base field.

        EXAMPLES:
            sage: k.<a> = NumberField([x^2 + 3, x^2 + 1])
            sage: m = k.base_field(); m
            Number Field in a1 with defining polynomial x^2 + 1
            sage: k._coerce_(m.0 + 2/3)
            a1 + 2/3
            sage: s = k._coerce_(m.0); s
            a1
            sage: s^2
            -1

        This implicitly tests this coercion map:
            sage: K.<a> = NumberField([x^2 + p for p in [5,3,2]])
            sage: K._coerce_(K.base_field().0)
            a1
            sage: K._coerce_(K.base_field().0)^2
            -3
        """
        abs_base, from_abs_base, to_abs_base = self.absolute_base_field()
        # Write element in terms of the absolute base field
        element = self.base_field()._coerce_impl(element)
        element = to_abs_base(element)
        # Obtain the polynomial in y corresponding to element in terms of the absolute base
        f = element.polynomial('y')
        # Find an expression in terms of the absolute generator for self of element.
        expr_x = self.pari_rnf().rnfeltreltoabs(f._pari_())
        # Convert to a SAGE polynomial, then to one in gen(), and return it
        R = self.polynomial_ring()
        return self(R(expr_x))

    def _fractional_ideal_class_(self):
        """
        Return the Python class used to represent ideals of a relative
        number field.

        EXAMPLES:
            sage: k.<a> = NumberField([x^5 + 2, x^7 + 3])
            sage: k._fractional_ideal_class_ ()
            <class 'sage.rings.number_field.number_field_ideal_rel.NumberFieldFractionalIdeal_rel'>
        """
        return sage.rings.number_field.number_field_ideal_rel.NumberFieldFractionalIdeal_rel

    def _pari_base_bnf(self, proof=None):
        """
        Return the PARI bnf (big number field) representation of the
        base field.

        INPUT:
            proof -- bool (default: True) if True, certify correctness
                     of calculations (not assuming GRH).

        EXAMPLES:
            sage: k.<a> = NumberField([x^3 + 2, x^2 + 2])
            sage: k._pari_base_bnf()
            [[;], matrix(0,9), [;], ... 0]
        """
        proof = proof_flag(proof)
        # No need to certify the same field twice, so we'll just check
        # that the base field is certified.
        if proof:
            self.base_field().pari_bnf_certify()
        return self.__base_bnf

    def _pari_base_nf(self):
        """
        Return the PARI number field representation of the base field.

        EXAMPLES:
            sage: y = polygen(QQ,'y')
            sage: k.<a> = NumberField([y^3 + 2, y^2 + 2])
            sage: k._pari_base_nf()
            [y^2 + 2, [0, 1], -8, 1, ..., [1, 0, 0, -2; 0, 1, 1, 0]]
        """
        return self.__base_nf

    def is_galois(self):
        r"""
        Return True if this relative number field is Galois over $\QQ$.

        EXAMPLES:
            sage: k.<a> =NumberField([x^3 - 2, x^2 + x + 1])
            sage: k.is_galois()
            True
            sage: k.<a> =NumberField([x^3 - 2, x^2 + 1])
            sage: k.is_galois()
            False
        """
        return self.absolute_field('a').is_galois()

    def vector_space(self):
        """
        Return vector space over the base field of self and isomorphisms
        from the vector space to self and in the other direction.

        EXAMPLES:
            sage: K.<a,b,c> = NumberField([x^2 + 2, x^3 + 2, x^3 + 3]); K
            Number Field in a with defining polynomial x^2 + 2 over its base field
            sage: V, from_V, to_V = K.vector_space()
            sage: from_V(V.0)
            1
            sage: to_V(K.0)
            (0, 1)
            sage: from_V(to_V(K.0))
            a
            sage: to_V(from_V(V.0))
            (1, 0)
            sage: to_V(from_V(V.1))
            (0, 1)

        The underlying vector space and maps is cached:
            sage: W, from_V, to_V = K.vector_space()
            sage: V is W
            True
        """
        try:
            return self.__vector_space
        except AttributeError:
            pass
        V = self.base_field()**self.degree()
        from_V = maps.MapRelativeVectorSpaceToRelativeNumberField(V, self)
        to_V   = maps.MapRelativeNumberFieldToRelativeVectorSpace(self, V)
        self.__vector_space = (V, from_V, to_V)
        return self.__vector_space

    def absolute_vector_space(self):
        """
        EXAMPLES:
            sage: K.<a,b> = NumberField([x^3 + 3, x^3 + 2]); K
            Number Field in a with defining polynomial x^3 + 3 over its base field
            sage: V,from_V,to_V = K.absolute_vector_space(); V
            Vector space of dimension 9 over Rational Field
            sage: from_V
            Isomorphism from Vector space of dimension 9 over Rational Field to Number Field in a with defining polynomial x^3 + 3 over its base field
            sage: to_V
            Isomorphism from Number Field in a with defining polynomial x^3 + 3 over its base field to Vector space of dimension 9 over Rational Field
            sage: c = (a+1)^5; c
            7*a^2 - 10*a - 29
            sage: to_V(c)
            (-29, -712/9, 19712/45, 0, -14/9, 364/45, 0, -4/9, 119/45)
            sage: from_V(to_V(c))
            7*a^2 - 10*a - 29
            sage: from_V(3*to_V(b))
            3*b
        """
        try:
            return self.__absolute_vector_space
        except AttributeError:
            pass
        K = self.absolute_field('a')
        from_K, to_K = K.structure()
        V, from_V, to_V = K.vector_space()
        fr = maps.MapVectorSpaceToRelativeNumberField(V, self, from_V, from_K)
        to   = maps.MapRelativeNumberFieldToVectorSpace(self, V, to_K, to_V)
        ans = (V, fr, to)
        self.__absolute_vector_space = ans
        return ans

    def absolute_base_field(self):
        """
        Return the base field of this relative extension, but viewed
        as an absolute field over QQ.

        EXAMPLES:
            sage: K.<a,b,c> = NumberField([x^2 + 2, x^3 + 3, x^3 + 2])
            sage: K
            Number Field in a with defining polynomial x^2 + 2 over its base field
            sage: K.base_field()
            Number Field in b with defining polynomial x^3 + 3 over its base field
            sage: K.absolute_base_field()[0]
            Number Field in a with defining polynomial x^9 + 3*x^6 + 165*x^3 + 1
            sage: K.base_field().absolute_field('z')
            Number Field in z with defining polynomial x^9 + 3*x^6 + 165*x^3 + 1
        """
        return self.__absolute_base_field

    def _gen_relative(self):
        """
        Return root of defining polynomial, which is a generator of
        the relative number field over the base.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2+1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: y = polygen(k)
            sage: m.<b> = k.extension(y^2+3); m
            Number Field in b with defining polynomial x^2 + 3 over its base field
            sage: c = m.gen(); c
            b
            sage: c^2 + 3
            0
        """
        try:
            return self.__gen_relative
        except AttributeError:
            rnf = self.pari_rnf()
            f = (pari('x') - rnf[10][2]*rnf[10][1]).lift()
            self.__gen_relative = self._element_class(self, f)
            return self.__gen_relative

    def pari_polynomial(self):
        """
        PARI polynomial corresponding to the polynomial over the
        rationals that defines this field as an absolute number field.

        EXAMPLES:
            sage: k.<a, c> = NumberField([x^2 + 3, x^2 + 1])
            sage: k.pari_polynomial()
            x^4 + 8*x^2 + 4
            sage: k.defining_polynomial ()
            x^2 + 3
        """
        try:
            return self.__pari_polynomial
        except AttributeError:
            poly = self.absolute_polynomial()
            with localvars(poly.parent(), 'x'):
                self.__pari_polynomial = poly._pari_()
            return self.__pari_polynomial

    def pari_rnf(self):
        """
        Return the PARI relative number field object associated
        to this relative extension.

        EXAMPLES:
            sage: k.<a> = NumberField([x^4 + 3, x^2 + 2])
            sage: k.pari_rnf()
            [x^4 + 3, [], [[108, 0; 0, 108], [3, 0]~], ... 0]
        """
        return self.__rnf

    def pari_relative_polynomial(self):
        """
        Return the PARI relative polynomial associated to this
        number field.  This is always a polynomial in x and y.

        EXAMPLES:
            sage: k.<i> = NumberField(x^2 + 1)
            sage: m.<z> = k.extension(k['w']([i,0,1]))
            sage: m
            Number Field in z with defining polynomial w^2 + i over its base field
            sage: m.pari_relative_polynomial ()
            x^2 + y
        """
        return self.__pari_relative_polynomial

    def number_of_roots_of_unity(self):
        """
        Return number of roots of unity in this relative field.

        EXAMPLES:
            sage: K.<a, b> = NumberField( [x^2 + x + 1, x^4 + 1] )
            sage: K.number_of_roots_of_unity()
            24
            sage: K.roots_of_unity()[:5]
            [-b^3*a, b^2*a + b^2, -b, -a, -b^3*a - b^3]
        """
        return self.absolute_field('a').number_of_roots_of_unity()

    def roots_of_unity(self):
        """
        Return all the roots of unity in this relative field, primitive or not.

        EXAMPLES:
            sage: K.<a, b> = NumberField( [x^2 + x + 1, x^4 + 1] )
            sage: K.roots_of_unity()[:5]
            [-b^3*a, b^2*a + b^2, -b, -a, -b^3*a - b^3]
        """
        abs = self.absolute_field('a')
        from_abs, _ = abs.structure()
        return [from_abs(x) for x in abs.roots_of_unity()]

    def absolute_generator(self):
        """
        Return the chosen generator over QQ for this relative number field.

        EXAMPLES:
            sage: y = polygen(QQ,'y')
            sage: k.<a> = NumberField([y^2 + 2, y^4 + 3])
            sage: g = k.absolute_generator(); g
            a0 - a1
            sage: g.minpoly()
            x^2 + 2*a1*x + a1^2 + 2
            sage: g.absolute_minpoly()
            x^8 + 8*x^6 + 30*x^4 - 40*x^2 + 49
        """
        try:
            return self.__abs_gen
        except AttributeError:
            self.__abs_gen = self._element_class(self, QQ['x'].gen())
            return self.__abs_gen


    def absolute_field(self, names):
        r"""
        Return an absolute number field K that is isomorphic to this
        field along with a field-theoretic bijection from self to K
        and from K to self.

        INPUT:
            names -- string; name of generator of the absolute field

        OUTPUT:
            K -- an absolute number field

        Also, \code{K.structure()} returns from_K and to_K, where
        from_K is an isomorphism from K to self and to_K is an isomorphism
        from self to K.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: L.<xyz> = K.absolute_field(); L
            Number Field in xyz with defining polynomial x^8 + 8*x^6 + 30*x^4 - 40*x^2 + 49
            sage: L.<c> = K.absolute_field(); L
            Number Field in c with defining polynomial x^8 + 8*x^6 + 30*x^4 - 40*x^2 + 49

            sage: from_L, to_L = L.structure()
            sage: from_L
            Isomorphism from Number Field in c with defining polynomial x^8 + 8*x^6 + 30*x^4 - 40*x^2 + 49 to Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: from_L(c)
            a - b
            sage: to_L
            Isomorphism from Number Field in a with defining polynomial x^4 + 3 over its base field to Number Field in c with defining polynomial x^8 + 8*x^6 + 30*x^4 - 40*x^2 + 49
            sage: to_L(a)
            -5/182*c^7 - 87/364*c^5 - 185/182*c^3 + 323/364*c
            sage: to_L(b)
            -5/182*c^7 - 87/364*c^5 - 185/182*c^3 - 41/364*c
            sage: to_L(a)^4
            -3
            sage: to_L(b)^2
            -2
        """
        try:
            return self.__absolute_field[names]
        except KeyError:
            pass
        except AttributeError:
            self.__absolute_field = {}
        K = NumberField(self.absolute_polynomial(), names, cache=False)
        from_K = maps.MapAbsoluteToRelativeNumberField(K, self)
        to_K = maps.MapRelativeToAbsoluteNumberField(self, K)
        K._set_structure(from_K, to_K)
        self.__absolute_field[names] = K
        return K

    def absolute_polynomial_ntl(self):
        """
        Return defining polynomial of this number field
        as a pair, an ntl polynomial and a denominator.

        This is used mainly to implement some internal arithmetic.

        EXAMPLES:
            sage: NumberField(x^2 + (2/3)*x - 9/17,'a').polynomial_ntl()
            ([-27 34 51], 51)
        """
        try:
            return (self.__abs_polynomial_ntl, self.__abs_denominator_ntl)
        except AttributeError:
            self.__abs_denominator_ntl = ntl.ZZ()
            den = self.absolute_polynomial().denominator()
            self.__abs_denominator_ntl.set_from_sage_int(ZZ(den))
            self.__abs_polynomial_ntl = ntl.ZZX((self.absolute_polynomial()*den).list())
        return (self.__abs_polynomial_ntl, self.__abs_denominator_ntl)

    def absolute_polynomial(self):
        r"""
        Return the polynomial over $\QQ$ that defines this field as an
        extension of the rational numbers.

        EXAMPLES:
            sage: k.<a, b> = NumberField([x^2 + 1, x^3 + x + 1]); k
            Number Field in a with defining polynomial x^2 + 1 over its base field
            sage: k.absolute_polynomial()
            x^6 + 5*x^4 - 2*x^3 + 4*x^2 + 4*x + 1
        """
        try:
            return self.__absolute_polynomial
        except AttributeError:
            pbn = self._pari_base_nf()
            prp = self.pari_relative_polynomial()
            pari_poly = pbn.rnfequation(prp)
            R = QQ['x']
            self.__absolute_polynomial = R(pari_poly)
            return self.__absolute_polynomial

    def base_field(self):
        """
        Return the base field of this relative number field.

        EXAMPLES:
            sage: k.<a> = NumberField([x^3 + x + 1])
            sage: R.<z> = k[]
            sage: L.<b> = NumberField(z^3 + a)
            sage: L.base_field()
            Number Field in a with defining polynomial x^3 + x + 1
            sage: L.base_field() is k
            True

        This is very useful because the print representation of
        a relative field doesn't describe the base field.
            sage: L
            Number Field in b with defining polynomial z^3 + a over its base field
        """
        return self.__base_field

    def base_ring(self):
        """
        This is exactly the same as base_field.

        EXAMPLES:
            sage: k.<a> = NumberField([x^2 + 1, x^3 + x + 1])
            sage: k.base_ring()
            Number Field in a1 with defining polynomial x^3 + x + 1
            sage: k.base_field()
            Number Field in a1 with defining polynomial x^3 + x + 1
        """
        return self.base_field()

    def embeddings(self, K):
        """
        Compute all field embeddings of the relative number field self
        into the field K (which need not even be a number field, e.g.,
        it could be the complex numbers). This will return an
        identical result when given K as input again.

        If possible, the most natural embedding of K into self
        is put first in the list.

        INPUT:
            K -- a number field

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^3 - 2, x^2+1])
            sage: f = K.embeddings(ComplexField(58)); f
            [
            Relative number field morphism:
              From: Number Field in a with defining polynomial x^3 - 2 over its base field
              To:   Complex Field with 58 bits of precision
              Defn: a |--> -0.62996052494743676 - 1.0911236359717214*I
                    b |--> -1.9428902930940239e-16 + 1.0000000000000000*I,
            ...
            Relative number field morphism:
              From: Number Field in a with defining polynomial x^3 - 2 over its base field
              To:   Complex Field with 58 bits of precision
              Defn: a |--> 1.2599210498948731
                    b |--> -0.99999999999999999*I
            ]
            sage: f[0](a)^3
            2.0000000000000002 - 8.6389229103644993e-16*I
            sage: f[0](b)^2
            -1.0000000000000001 - 3.8857805861880480e-16*I
            sage: f[0](a+b)
            -0.62996052494743693 - 0.091123635971721295*I
        """
        try:
            return self.__embeddings[K]
        except AttributeError:
            self.__embeddings = {}
        except KeyError:
            pass
        L = self.absolute_field('a')
        E = L.embeddings(K)
        v = [self.hom(f, K) for f in E]

        # If there is an embedding that preserves variable names
        # then it is most natural, so we put it first.
        put_natural_embedding_first(v)

        self.__embeddings[K] = Sequence(v, cr=v!=[], immutable=True, check=False, universe=self.Hom(K))
        return self.__embeddings[K]

    def relative_discriminant(self, proof=None):
        r"""
        Return the relative discriminant of this extension $L/K$ as
        an ideal of $K$.  If you want the (rational) discriminant of
        $L/Q$, use e.g. \code{L.discriminant()}.

        TODO: Note that this uses PARI's \code{rnfdisc} function, which
        according to the documentation takes an \code{nf} parameter in
        GP but a \code{bnf} parameter in the C library.  If the C
        library actually accepts an \code{nf}, then this function
        should be fixed and the \code{proof} parameter removed.

        INPUT:
            proof -- (default: False)

        EXAMPLE:
            sage: K.<i> = NumberField(x^2 + 1)
            sage: t = K['t'].gen()
            sage: L.<b> = K.extension(t^4 - i)
            sage: L.relative_discriminant()
            Fractional ideal (256)
            sage: factor(L.discriminant())
            2^24
            sage: factor( L.relative_discriminant().norm() )
            2^16
        """
        proof = proof_flag(proof)

        bnf = self._pari_base_bnf(proof)
        K = self.base_field()
        R = K.polynomial().parent()
        D, d = bnf.rnfdisc(self.pari_relative_polynomial())
        return K.ideal([ K(R(x)) for x in convert_from_zk_basis(K, D) ])

    def order(self, *gens, **kwds):
        """
        Return the order with given ring generators in the maximal
        order of this number field.

        INPUT:
            gens -- list of elements of self; if no generators are
                    given, just returns the cardinality of this number
                    field (oo) for consistency.
            check_is_integral -- bool (default: True), whether to check
                  that each generator is integral.
            check_rank -- bool (default: True), whether to check that
                  the ring generated by gens is of full rank.
            allow_subfield -- bool (default: False), if True and the generators
                  do not generate an order, i.e., they generate a subring
                  of smaller rank, instead of raising an error, return
                  an order in a smaller number field.

        The check_is_integral and check_rank inputs must be given as
        explicit keyword arguments.

        EXAMPLES:
            sage: P.<a,b,c> = QQ[2^(1/2), 2^(1/3), 3^(1/2)]
            sage: R = P.order([a,b,c]); R
            Relative Order in Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field

        The base ring of an order in a relative extension is still ZZ.
            sage: R.base_ring()
            Integer Ring

        One must give enough generators to generate a ring of finite index
        in the maximal order:
            sage: P.order([a,b])
            Traceback (most recent call last):
            ...
            ValueError: the rank of the span of gens is wrong
        """
        import sage.rings.number_field.order as order
        if len(gens) == 0:
            return NumberField_generic.order(self)
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        gens = [self(x) for x in gens]
        return order.relative_order_from_ring_generators(gens, **kwds)


    def galois_group(self, pari_group = True, algorithm='pari'):
        r"""
        Return the Galois group of the Galois closure of this number
        field as an abstract group.  Note that even though this is an
        extension $L/K$, the group will be computed as if it were $L/\QQ$.

        INPUT:
            pari_group -- bool (default: False); if True instead return
                          the Galois group as a PARI group.
            algorithm -- 'pari', 'kash', 'magma' (default: 'pari', except
                          when the degree is >= 12 when 'kash' is tried)

        For more (important!) documentation, so the documentation
        for Galois groups of polynomials over $\QQ$, e.g., by
        typing \code{K.polynomial().galois_group?}, where $K$
        is a number field.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^2 + 1)
            sage: R.<t> = PolynomialRing(K)
            sage: L = K.extension(t^5-t+a, 'b')
            sage: L.galois_group()
            Galois group PARI group [240, -1, 22, "S(5)[x]2"] of degree 10 of the Number Field in b with defining polynomial t^5 - t + a over its base field
        """
        try:
            return self.__galois_group[pari_group, algorithm]
        except KeyError:
            pass
        except AttributeError:
            self.__galois_group = {}

        G = self.absolute_polynomial().galois_group(pari_group = pari_group,
                                                    algorithm = algorithm)
        H = GaloisGroup(G, self)
        self.__galois_group[pari_group, algorithm] = H
        return H


    def is_free(self, proof=None):
        r"""
        Determine whether or not $L/K$ is free (i.e. if $\mathcal{O}_L$ is
        a free $\mathcal{O}_K$-module).

        INPUT:
            proof -- default: True

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^2+6)
            sage: L.<b> = K.extension(K['x'].gen()^2 + 3)    ## extend by x^2+3
            sage: L.is_free()
            False
        """
        proof = proof_flag(proof)
        base_bnf = self._pari_base_bnf(proof)
        if base_bnf.rnfisfree(self.pari_relative_polynomial()) == 1:
            return True
        return False

    def lift_to_base(self, element):
        """
        Lift an element of this extension into the base field if possible,
        or raise a ValueError if it is not possible.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^3 - 2, 'a')
            sage: R = K['x']
            sage: L = K.extension(R.gen()^2 - K.gen(), 'b')
            sage: b = L.gen()
            sage: L.lift_to_base(b^4)
            a^2
            sage: L.lift_to_base(b)
            Traceback (most recent call last):
            ...
            ValueError: The element b is not in the base field
        """
        poly_xy = self.pari_rnf().rnfeltabstorel( self(element)._pari_() )
        str_poly = str(poly_xy)
        if str_poly.find('x') >= 0:
            raise ValueError, "The element %s is not in the base field"%element
        # We convert to a string to avoid some serious nastiness with
        # PARI polynomials secretely thinkining they are in more variables
        # than they are.
        f = QQ['y'](str_poly)
        return self.base_field()(f.list())

    def polynomial(self):
        """
        Return the defining polynomial of this number field.

        EXAMPLES:
            sage: y = polygen(QQ,'y')
            sage: k.<a> = NumberField([y^2 + y + 1, x^3 + x + 1])
            sage: k.polynomial()
            y^2 + y + 1

        This is the same as defining_polynomial:
            sage: k.defining_polynomial()
            y^2 + y + 1

        Use absolute polynomial for a polynomial that defines the
        absolute extension.
            sage: k.absolute_polynomial()
            x^6 + 3*x^5 + 8*x^4 + 9*x^3 + 7*x^2 + 6*x + 3
        """
        return self.__relative_polynomial

    def relativize(self, alpha, names):
        r"""
        Given an element alpha in self, return a relative number field
        $K$ isomorphic to self that is relative over the absolute field
        $\QQ(\alpha)$, along with isomorphisms from $K$ to self and
        from self to K.

        INPUT:
            alpha -- an element of self.
            names -- name of generator for output field K.

        OUTPUT:
            K -- relative number field

        Also, \code{K.structure()} returns from_K and to_K, where
        from_K is an isomorphism from K to self and to_K is an isomorphism
        from self to K.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^4 + 3, x^2 + 2]); K
            Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: L.<z,w> = K.relativize(a^2)
            sage: z^2
            z^2
            sage: w^2
            -3
            sage: L
            Number Field in z with defining polynomial x^4 + (-2*w + 4)*x^2 + 4*w + 1 over its base field
            sage: L.base_field()
            Number Field in w with defining polynomial x^2 + 3
        """
        K = self.absolute_field('a')
        from_K, to_K = K.structure()
        beta = to_K(alpha)
        S = K.relativize(beta, names)
        # Now S is the appropriate field,
        # but the structure maps attached to S
        # are isomorphisms with the absolute
        # field.  We have to compose them
        # with from_K and to_K to get
        # the appropriate maps.
        from_S, to_S = S.structure()

        # Map from S to self:
        #   x |--> from_K(from_S(x))
        # Map from self to S:
        #   x |--> to_K(from_K(x))
        new_to_S = self.Hom(S)(to_S)
        a = from_S.abs_hom()
        W = a.domain()
        phi = W.hom([from_K(a(W.gen()))])
        new_from_S = S.Hom(self)(phi)
        S._set_structure(new_from_S, new_to_S, unsafe_force_change=True)
        return S

class NumberField_cyclotomic(NumberField_absolute):
    """
    Create a cyclotomic extension of the rational field.

    The command CyclotomicField(n) creates the n-th cyclotomic field,
    obtained by adjoining an n-th root of unity to the rational field.

    EXAMPLES:
        sage: CyclotomicField(3)
        Cyclotomic Field of order 3 and degree 2
        sage: CyclotomicField(18)
        Cyclotomic Field of order 18 and degree 6
        sage: z = CyclotomicField(6).gen(); z
        zeta6
        sage: z^3
        -1
        sage: (1+z)^3
        6*zeta6 - 3

        sage: K = CyclotomicField(197)
        sage: loads(K.dumps()) == K
        True
        sage: loads((z^2).dumps()) == z^2
        True

        sage: cf12 = CyclotomicField( 12 )
        sage: z12 = cf12.0
        sage: cf6 = CyclotomicField( 6 )
        sage: z6 = cf6.0
        sage: FF = Frac( cf12['x'] )
        sage: x = FF.0
        sage: z6*x^3/(z6 + x)
        zeta12^2*x^3/(x + zeta12^2)

        sage: cf6 = CyclotomicField(6) ; z6 = cf6.gen(0)
        sage: cf3 = CyclotomicField(3) ; z3 = cf3.gen(0)
        sage: cf3(z6)
        zeta3 + 1
        sage: cf6(z3)
        zeta6 - 1
        sage: type(cf6(z3))
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
        sage: cf1 = CyclotomicField(1) ; z1 = cf1.0
        sage: cf3(z1)
        1
        sage: type(cf3(z1))
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
    """
    def __init__(self, n, names):
        """
        A cyclomotic field, i.e., a field obtained by adjoining an
        n-th root of unity to the rational numbers.

        EXAMPLES:
            sage: k = CyclotomicField(3)
            sage: type(k)
            <class 'sage.rings.number_field.number_field.NumberField_cyclotomic'>
        """
        f = QQ['x'].cyclotomic_polynomial(n)
        if names[0].startswith('zeta'):
            latex_name = "\\zeta_{%s}"%n
        else:
            latex_name = None
        NumberField_absolute.__init__(self, f,
                                      name= names,
                                      latex_name=latex_name,
                                      check=False)
        n = integer.Integer(n)
        self.__n = n
        if n%2:
            self.__zeta_order = 2*n
        else:
            self.__zeta_order = n
        ## quadratic number fields require this:
        if f.degree() == 2:
            self._element_class = number_field_element_quadratic.NumberFieldElement_quadratic
            if n == 4:
                self._D = ZZ(-1)
                self._NumberField_generic__gen = self._element_class(self, (QQ(0), QQ(1)))
            else:
                ## n is 3 or 6
                self._D = ZZ(-3)
                one_half = ZZ(1)/ZZ(2)
                if n == 3:
                    self._NumberField_generic__gen = self._element_class(self, (one_half-1, one_half))
                else:
                    self._NumberField_generic__gen = self._element_class(self, (one_half, one_half))
        zeta = self.gen()
        zeta._set_multiplicative_order(n)

    def __reduce__(self):
        """
        TESTS:
            sage: K.<zeta7> = CyclotomicField(7)
            sage: L = loads(dumps(K))
            sage: L
            Cyclotomic Field of order 7 and degree 6
            sage: L == K
            True
        """
        return NumberField_cyclotomic_v1, (self.__n, self.variable_name())

    def _magma_init_(self):
        # TODO: I really don't like this on multiple levels.
        # (1) it kills a global symbol self.gen()
        # (2) it abuses how conversion works and throws in an extra define.
        # (3) a cyclo field in a funny generator wouldn't get converted to
        #     one with the right name via this.
        # (4) One should define _magma_convert_ instead of _magma_init_
        #     in this case, probably.
        #    -- William
        """
        Function returning a string to create this cyclotomic field in Magma.

        NOTE: The Magma generator name is also initialized to be the
        same as for the Sage field.

        EXAMPLES:
            sage: K=CyclotomicField(7,'z')
            sage: K._magma_init_()                                # optional - magma
            'CyclotomicField(7); z:=CyclotomicField(7).1;'
            sage: K=CyclotomicField(7,'zeta')
            sage: K._magma_init_()                                # optional - magma
            'CyclotomicField(7); zeta:=CyclotomicField(7).1;'
        """
        return 'CyclotomicField(%s); %s:=CyclotomicField(%s).1;'%(self.__n, self.gen(), self.__n)

    def _repr_(self):
        r"""
        Return string representation of this cyclotomic field.

        The ``order'' of the cyclotomic field $\QQ(\zeta_n)$ in the
        string output refers to the order of the $\zeta_n$, i.e., it
        is the integer $n$.  The degree is the degree of the field as
        an extension of $\QQ$.

        EXAMPLES:
            sage: CyclotomicField(4)._repr_()
            'Cyclotomic Field of order 4 and degree 2'
            sage: CyclotomicField(400)._repr_()
            'Cyclotomic Field of order 400 and degree 160'
        """
        return "Cyclotomic Field of order %s and degree %s"%(
                self.__n, self.degree())

    def _n(self):
        """
        Return the n used to create this cyclotomic field.

        EXAMPLES:
            sage: CyclotomicField(3).zeta_order()
            6
            sage: CyclotomicField(3)._n()
            3
        """
        return self.__n

    def _latex_(self):
        """
        Return the latex representation of this cyclotomic field.

        EXAMPLES:
            sage: Z = CyclotomicField(4)
            sage: Z.gen()
            zeta4
            sage: latex(Z)
            \mathbf{Q}(\zeta_{4})

        Latex printing respects the generator name.
            sage: k.<a> = CyclotomicField(4)
            sage: latex(k)
            \mathbf{Q}[a]/(a^{2} + 1)
            sage: k
            Cyclotomic Field of order 4 and degree 2
            sage: k.gen()
            a
        """
        v = self.latex_variable_name()
        if v.startswith('\\zeta_'):
            return "%s(%s)"%(latex(QQ), v)
        else:
            return NumberField_generic._latex_(self)

    def __call__(self, x):
        """
        Create an element of this cyclotomic field from $x$.

        EXAMPLES:
        The following example illustrates coercion from the cyclotomic
        field Q(zeta_42) to the cyclotomic field Q(zeta_6), in a case
        where such coercion is defined:

            sage: k42 = CyclotomicField(42)
            sage: k6 = CyclotomicField(6)
            sage: a = k42.gen(0)
            sage: b = a^7
            sage: b
            zeta42^7
            sage: k6(b)
            zeta6
            sage: b^2
            zeta42^7 - 1
            sage: k6(b^2)
            zeta6 - 1

        Coercion of GAP cyclotomic elements is also supported.

        EXAMPLE:
            sage: K.<z> = CyclotomicField(7)
            sage: O = K.maximal_order()
            sage: K(O.1)
            z
            sage: K(O.1^2 + O.1 - 2)
            z^2 + z - 2
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            if isinstance(x.parent(), NumberField_cyclotomic):
                return self._coerce_from_other_cyclotomic_field(x)
            else:
                return NumberField_absolute.__call__(self, x)
        elif sage.interfaces.gap.is_GapElement(x):
            return self._coerce_from_gap(x)
        elif isinstance(x,str):
            return self._coerce_from_str(x)
        else:
            return self._coerce_non_number_field_element_in(x)

    # TODO:
    # The following is very nice and much more flexible / powerful.
    # However, it is simply not *consistent*, since it totally
    # breaks the doctests in eisenstein_submodule.py.
    # FIX THIS.

##     def _will_be_better_coerce_from_other_cyclotomic_field(self, x, only_canonical=False):
##         """
##         Coerce an element x of a cyclotomic field into self, if at all possible.

##         INPUT:
##             x -- number field element

##             only_canonical -- bool (default: False); Attempt to work,
##                    even in some cases when x is not in a subfield of
##                    the cyclotomics (as long as x is a root of unity).

##         EXAMPLES:
##             sage: k5 = CyclotomicField(5)
##             sage: k3 = CyclotomicField(3)
##             sage: k15 = CyclotomicField(15)
##             sage: k15._coerce_from_other_cyclotomic_field(k3.gen())
##             zeta15^5
##             sage: k15._coerce_from_other_cyclotomic_field(k3.gen()^2 + 17/3)
##             -zeta15^5 + 14/3
##             sage: k3._coerce_from_other_cyclotomic_field(k15.gen()^5)
##             zeta3
##             sage: k3._coerce_from_other_cyclotomic_field(-2/3 * k15.gen()^5 + 2/3)
##             -2/3*zeta3 + 2/3
##         """

##         K = x.parent()

##         if K is self:
##             return x
##         elif K == self:
##             return self._element_class(self, x.polynomial())
##         n = K.zeta_order()
##         m = self.zeta_order()
##         print n, m, x


##         self_gen = self.gen()

##         if m % n == 0:   # easy case
##             # pass this off to a method in the element class
##             # it can be done very quickly and easily by the cython<->NTL
##             # interface there
##             return x._lift_cyclotomic_element(self)

##         # Whatever happens below, it has to be consistent with
##         #  zeta_r |---> (zeta_s)^m

##         if m % 2 and not n%2:
##             m *= 2
##             self_gen = -self_gen

##         if only_canonical and m % n:
##             raise TypeError, "no canonical coercion"

##         if not is_CyclotomicField(K):
##             raise TypeError, "x must be in a cyclotomic field"

##         v = x.list()

##         # Find the smallest power r >= 1 of the generator g of K that is in self,
##         # i.e., find the smallest r such that g^r has order dividing m.

##         d = sage.rings.arith.gcd(m,n)
##         r = n // d

##         # Since we use the power basis for cyclomotic fields, if every
##         # v[i] with i not divisible by r is 0, then we're good.

##         # If h generates self and has order m, then the element g^r
##         # maps to the power of self of order gcd(m,n)., i.e., h^(m/gcd(m,n))
##         #
##         z = self_gen**(m // d)
##         w = self(1)

##         a = self(0)
##         for i in range(len(v)):
##             if i%r:
##                 if v[i]:
##                     raise TypeError, "element does not belong to cyclotomic field"
##             else:
##                 a += w*v[i]
##                 w *= z
##         return a

    def _coerce_from_other_cyclotomic_field(self, x, only_canonical=False):
        """
        Coerce an element x of a cyclotomic field into self, if at all possible.

        INPUT:
            x -- number field element
            only_canonical -- bool (default: False); Attempt to work, even in some
                   cases when x is not in a subfield of the cyclotomics (as long as x is
                   a root of unity).

        EXAMPLES:
            sage: K = CyclotomicField(24) ; L = CyclotomicField(48)
            sage: L._coerce_from_other_cyclotomic_field(K.0+1)
            zeta48^2 + 1
            sage: K(L.0**2)
            zeta24
        """
        K = x.parent()
        if K is self:
            return x
        elif K == self:
            return self._element_class(self, x.polynomial())
        n = K._n()
        m = self._n()
        if m % n == 0:   # easy case
            # pass this off to a method in the element class
            # it can be done very quickly and easily by the
            # Cython<->NTL interface there
            return x._lift_cyclotomic_element(self)
        else:
            if only_canonical:
                raise TypeError
            n = x.multiplicative_order()
            m = self.zeta_order()
            if m % n == 0:
                # Harder case.  E.g., x = (zeta_42)^7 and
                # self.__zeta = zeta_6, so it is possible to
                # coerce x in, but not zeta_42 in.
                # Algorithm:
                #    1. Compute self.__zeta as an element
                #       of K = parent of x.  Call this y.
                #    2. Write x as a power r of y.
                #       TODO: we do step two STUPIDLY.
                #    3. Return self.__zeta to the power r.
                y = K(self.zeta(m))
                z = y
                for r in xrange(y.multiplicative_order()):
                    if z == x:
                        return self.zeta(m)**(r+1)
                    z *= y
            raise TypeError, "Cannot coerce %s into %s"%(x,self)
        return self._element_class(self, g)


    def _coerce_from_gap(self, x):
        """
        Attempt to coerce a GAP number field element into this cyclotomic field.

        EXAMPLES:
            sage: k5.<z> = CyclotomicField(5)
            sage: gap('E(5)^7 + 3')
            -3*E(5)-2*E(5)^2-3*E(5)^3-3*E(5)^4
            sage: w = gap('E(5)^7 + 3')
            sage: z^7 + 3
            z^2 + 3
            sage: k5(w)
            z^2 + 3
        """
        s = str(x)
        i = s.find('E(')
        if i == -1:
            return self(rational.Rational(s))
        j = i + s[i:].find(')')
        n = int(s[i+2:j])
        if n == self.zeta_order():
            K = self
        else:
            K = CyclotomicField(n)
        zeta = K.gen()
        s = s.replace('E(%s)'%n,'zeta')
        s = sage.misc.all.sage_eval(s, locals={'zeta':K.gen()})
        if K is self:
            return s
        else:
            return self(s)

    def _coerce_impl(self, x):
        """
        Canonical implicit coercion of x into self.

        Elements of other compatible cyclotomic fields coerce in, as
        do elements of the rings that coerce to all number fields
        (e.g., integers, rationals).

        EXAMPLES:
            sage: CyclotomicField(15)._coerce_impl(CyclotomicField(5).0 - 17/3)
            zeta15^3 - 17/3
            sage: K.<a> = CyclotomicField(16)
            sage: K(CyclotomicField(4).0)
            a^4
        """
        if isinstance(x, number_field_element.NumberFieldElement) and \
                isinstance(x.parent(), NumberField_cyclotomic):
            return self._coerce_from_other_cyclotomic_field(x, only_canonical=True)
        return NumberField_generic._coerce_impl(self, x)

    def _Hom_(self, codomain, cat=None):
        """
        Return homset of homomorphisms from the cyclotomic field self to the number field codomain.

        The cat option is currently ignored.

        EXAMPLES:
        This function is implicitly caled by the Hom method or function.
            sage: K.<a> = NumberField(x^2 + 3); K
            Number Field in a with defining polynomial x^2 + 3
            sage: CyclotomicField(3).Hom(K)
            Set of field embeddings from Cyclotomic Field of order 3 and degree 2 to Number Field in a with defining polynomial x^2 + 3
            sage: End(CyclotomicField(21))
            Automorphism group of Cyclotomic Field of order 21 and degree 12
        """
        import morphism
        return morphism.CyclotomicFieldHomset(self, codomain)

    def is_galois(self):
        """
        Return True since all cyclotomic fields are automatically Galois.

        EXAMPLES:
            sage: CyclotomicField(29).is_galois()
            True
        """
        return True

    def is_isomorphic(self, other):
       """
       Return True if the cyclotomic field self is isomorphic
       as a number field to other.

       EXAMPLES:
           sage: CyclotomicField(11).is_isomorphic(CyclotomicField(22))
           True
           sage: CyclotomicField(11).is_isomorphic(CyclotomicField(23))
           False
           sage: CyclotomicField(3).is_isomorphic(NumberField(x^2 + x +1, 'a'))
           True
       """
       if not isinstance(other, NumberField_generic):
           raise ValueError, "other must be a generic number field."
       return self._Hom_(other).order() > 0

    def complex_embedding(self, prec=53):
        r"""
        Return the embedding of this cyclotomic field into the
        approximate complex field with precision prec obtained by
        sending the generator $\zeta$ of self to exp(2*pi*i/n), where
        $n$ is the multiplicative order of $\zeta$.

        If prec is 53 (the default), then the complex double field is
        used; otherwise the arbitrary precision (but slow) complex
        field is used.
        EXAMPLES:
            sage: C = CyclotomicField(4)
            sage: C.complex_embedding()
            Ring morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Complex Double Field
              Defn: zeta4 |--> 6.12323399574e-17 + 1.0*I

        Note in the example above that the way zeta is computed (using
        sin and cosine in MPFR) means that only the prec bits of the
        number after the decimal point are valid.

            sage: K = CyclotomicField(3)
            sage: phi = K.complex_embedding(10)
            sage: phi(K.0)
            -0.50 + 0.87*I
            sage: phi(K.0^3)
            1.0
            sage: phi(K.0^3 - 1)
            0
            sage: phi(K.0^3 + 7)
            8.0
        """
        if prec == 53:
            CC = sage.rings.complex_double.CDF
        else:
            CC = sage.rings.complex_field.ComplexField(prec)
        return self.hom([CC.zeta(self._n())], check=False)

    def complex_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the
        approximate complex field with precision prec.

        If prec is 53 (the default), then the complex double field is
        used; otherwise the arbitrary precision (but slow) complex
        field is used.  If you want 53-bit arbitrary precision then
        do \code{self.embeddings(ComplexField(53))}.

        EXAMPLES:
            sage: CyclotomicField(5).complex_embeddings()
            [
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Double Field
              Defn: zeta5 |--> 0.309016994375 + 0.951056516295*I,
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Double Field
              Defn: zeta5 |--> -0.809016994375 + 0.587785252292*I,
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Double Field
              Defn: zeta5 |--> -0.809016994375 - 0.587785252292*I,
            Ring morphism:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Complex Double Field
              Defn: zeta5 |--> 0.309016994375 - 0.951056516295*I
            ]
        """
        if prec == 53:
            CC = sage.rings.complex_double.CDF
        else:
            CC = sage.rings.complex_field.ComplexField(prec)
        try:
            return self.__embeddings[CC]
        except AttributeError:
            self.__embeddings = {}
        except KeyError:
            pass
        n = self._n()
        z = CC.zeta(n)
        X = [m for m in range(n) if sage.rings.arith.gcd(m,n) == 1]
        v = [self.hom([z**n], check=False) for n in X]
        self.__embeddings[CC] = Sequence(v, cr=True, immutable=True,
                                         check=False, universe=self.Hom(CC))
        return self.__embeddings[CC]

    def real_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the
        approximate real field with precision prec.

        If prec is 53 (the default), then the real double field is
        used; otherwise the arbitrary precision (but slow) real field
        is used.

        Mostly, of course, there are no such embeddings.

        EXAMPLES:
            sage: CyclotomicField(4).real_embeddings()
            []
            sage: CyclotomicField(2).real_embeddings()
            [
            Ring morphism:
              From: Cyclotomic Field of order 2 and degree 1
              To:   Real Double Field
              Defn: -1 |--> -1.0
            ]
        """
        if prec == 53:
            K = sage.rings.real_double.RDF
        else:
            K = sage.rings.real_mpfr.RealField(prec)
        n = self._n()
        if n > 2:
            return Sequence([], cr=False, immutable=True,
                                         check=False, universe=self.Hom(K))
        else:
            return self.embeddings(K)

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this cyclotomic field, respectively.

        Trivial since, apart from QQ, cyclotomic fields are totally complex.

        EXAMPLES:
            sage: CyclotomicField(5).signature()
            (0, 2)
            sage: CyclotomicField(2).signature()
            (1, 0)
        """
        m = ZZ(self.degree())
        if m == 1:
            return (ZZ(1), ZZ(0))
        else:
            return (ZZ(0), ZZ(m/2))

    def discriminant(self, v=None):
        """
        Returns the discriminant of the ring of integers of the cyclotomic field self,
        or if v is specified, the determinant of the trace pairing
        on the elements of the list v.

        Uses the formula for the discriminant of a prime power cyclotomic field
        and Hilbert Theorem 88 on the discriminant of composita.

        INPUT:
            v (optional) -- list of element of this number field
        OUTPUT:
            Integer if v is omitted, and Rational otherwise.

        EXAMPLES:
            sage: CyclotomicField(20).discriminant()
            4000000
            sage: CyclotomicField(18).discriminant()
            -19683
        """
        if v == None:
            try:
                return self.__disc
            except AttributeError:
                n = self._n()
                deg = self.degree()
                d = ZZ(1) # so that CyclotomicField(1).disc() has the right type
                factors = n.factor()
                for f in factors:
                    p = f[0]
                    r = f[1]
                    e = (r*p - r - 1)*deg/(p-1)
                    d *= p**e
                sign = 1
                if len(factors) == 1 and (n == 4 or factors[0][0].mod(4) == 3):
                    sign = -1
                elif len(factors) == 2 and factors[0] == (2, 1) and factors[1][0].mod(4) == 3:
                    sign = -1
                self.__disc = sign*d
                return self.__disc
        else:
            return NumberField_generic.discriminant(self, v)


    def next_split_prime(self, p=2):
        """
        Return the next prime integer $p$ that splits completely in
        this cyclotomic field (and does not ramify).

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: K.next_split_prime(7)
            13
        """
        n = self._n()
        while True:
            p = sage.rings.arith.next_prime(p)
            if p % n == 1:
                return p

    def integral_basis(self, v=None):
        """
        Return a list of elements of this number field that are a basis
        for the full ring of integers.

        This field is cyclomotic, so this is a trivial computation,
        since the power basis on the generator is an integral basis.
        Thus the v parameter is ignored.

        EXAMPLES:
            sage: CyclotomicField(5).integral_basis()
            [1, zeta5, zeta5^2, zeta5^3]
        """
        try:
            return self._integral_basis_dict[tuple()]
        except KeyError:
            v = tuple()
            z = self.gen()
            a = self(1)
            B = []
            for n in xrange(self.degree()):
                B.append(a)
                a *= z
            self._integral_basis_dict[tuple()] = B
            return B

    def _compute_integral_basis(self, v=None):
        """
        Alias for self.integral_basis().

        EXAMPLES:
            sage: len(CyclotomicField(137)._compute_integral_basis())
            136
            sage: CyclotomicField(17).integral_basis() == CyclotomicField(17)._compute_integral_basis()
            True
        """
        return self.integral_basis()


    def zeta_order(self):
        """
        Return the order of the maximal root of unity contained in
        this cyclotomic field.

        EXAMPLES:
            sage: CyclotomicField(1).zeta_order()
            2
            sage: CyclotomicField(4).zeta_order()
            4
            sage: CyclotomicField(5).zeta_order()
            10
            sage: CyclotomicField(5)._n()
            5
            sage: CyclotomicField(389).zeta_order()
            778
        """
        return self.__zeta_order

    def _multiplicative_order_table(self):
        """
        Return a dictionary that maps powers of zeta to their order.
        This makes computing the orders of the elements of finite
        order in this field faster.

        EXAMPLES:
            sage: v = CyclotomicField(6)._multiplicative_order_table()
            sage: w = v.items(); w.sort(); w
            [(-1, 2), (1, 1), (-x, 3), (-x + 1, 6), (x - 1, 3), (x, 6)]
        """
        try:
            return self.__multiplicative_order_table
        except AttributeError:
            t = {}
            x = self(1)
            n = self.zeta_order()
            m = 0
            zeta = self.zeta(n)
            # todo: this desperately needs to be optimized!!!
            for i in range(n):
                t[x.polynomial()] = n//arith.GCD(m,n)   # multiplicative_order of (zeta_n)**m
                x *= zeta
                m += 1
            self.__multiplicative_order_table = t
            return t

    def zeta(self, n=None, all=False):
        """
        Returns an element of multiplicative order $n$ in this this
        cyclotomic field, if there is one.  Raises a ValueError if
        there is not.

        INPUT:
            n -- integer (default: None, returns element of maximal order)
            all -- bool (default: False) -- whether to return a list of
                        all n-th roots.

        OUTPUT:
            root of unity or list

        EXAMPLES:
            sage: k = CyclotomicField(7)
            sage: k.zeta()
            zeta7
            sage: k.zeta().multiplicative_order()
            7
            sage: k = CyclotomicField(49)
            sage: k.zeta().multiplicative_order()
            49
            sage: k.zeta(7).multiplicative_order()
            7
            sage: k.zeta()
            zeta49
            sage: k.zeta(7)
            zeta49^7

            sage: K.<a> = CyclotomicField(7)
            sage: K.zeta(14, all=True)
            [-a^4, -a^5, a^5 + a^4 + a^3 + a^2 + a + 1, -a, -a^2, -a^3]
            sage: K.<a> = CyclotomicField(10)
            sage: K.zeta(20, all=True)
            Traceback (most recent call last):
            ...
            ValueError: n (=20) does not divide order of generator

            sage: K.<a> = CyclotomicField(5)
            sage: K.zeta(4)
            Traceback (most recent call last):
            ...
            ValueError: n (=4) does not divide order of generator
            sage: v = K.zeta(5, all=True); v
            [a, a^2, a^3, -a^3 - a^2 - a - 1]
            sage: [b^5 for b in v]
            [1, 1, 1, 1]
        """
        if n is None:
            return self.gen()
        else:
            n = integer.Integer(n)
            z = self.gen()
            m = z.multiplicative_order()
            if n % 2 == 0 and m % 2 == 1:
                # In the n-th cyclotomic field, n odd, there are
                # actually 2*n-th roots of unity, so we include them.
                z = -z**((m+1)//2) # -z
                m = 2*m
            if m % n != 0:
                raise ValueError, "n (=%s) does not divide order of generator"%n
                # use generic method (factor cyclotomic polynomial)
                #  -- this is potentially really slow, so don't do it.
                #return field.Field.zeta(self, n, all=all)
            a = z**(m//n)
            if all:
                v = [a]
                b = a*a
                for i in range(2,n):
                    if sage.rings.arith.gcd(i, n) == 1:
                        v.append(b)
                    b = b * a
                return v
            else:
                return a

    def number_of_roots_of_unity(self):
        """
        Return number of roots of unity in this cyclotomic field.

        EXAMPLES:
            sage: K.<a> = CyclotomicField(21)
            sage: K.number_of_roots_of_unity()
            42
        """
        n = self._n()
        if n%2:
            n *= 2
        return n

    def roots_of_unity(self):
        """
        Return all the roots of unity in this cyclotomic field, primitive or not.

        EXAMPLES:
            sage: K.<a> = CyclotomicField(3)
            sage: zs = K.roots_of_unity(); zs
            [1, a, -a - 1, -1, -a, a + 1]
            sage: [ z**K.number_of_roots_of_unity() for z in zs ]
            [1, 1, 1, 1, 1, 1]
        """
        z = self.gen()
        n = self._n()
        v = [z**k for k in range(n)]
        if n%2:
            v += [-x for x in v]
        return v

class NumberField_quadratic(NumberField_absolute):
    """
    Create a quadratic extension of the rational field.

    The command QuadraticField(a) creates the field Q(sqrt(a)).

    EXAMPLES:
        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3
        sage: QuadraticField(-4, 'b')
        Number Field in b with defining polynomial x^2 + 4
    """
    def __init__(self, polynomial, name=None, check=True):
        """
        Create a quadratic number field.

        EXAMPLES:
            sage: k.<a> = QuadraticField(5, check=False); k
            Number Field in a with defining polynomial x^2 - 5

        Don't do this:
            sage: k.<a> = QuadraticField(4, check=False); k
            Number Field in a with defining polynomial x^2 - 4
        """
        NumberField_absolute.__init__(self, polynomial, name=name, check=check)
        self._element_class = number_field_element_quadratic.NumberFieldElement_quadratic
        c, b, a = [rational.Rational(t) for t in self.defining_polynomial().list()]
        # set the generator
        Dpoly = b*b - 4*a*c
        D = Dpoly.numer() * Dpoly.denom()
        # this could be done extremely quickly in cython
        for p in sage.rings.arith.primes(100):
            p2 = p*p
            while D % p2 == 0:
                D //= p2
        self._D = D
        parts = -b/(2*a), (Dpoly/D).sqrt()/(2*a)
        self._NumberField_generic__gen = self._element_class(self, parts)


    def coerce_map_from_impl(self, S):
        """
        EXAMPLES:
            sage: K.<a> = QuadraticField(-3)
            sage: f = K.coerce_map_from(QQ); f
            Natural morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 + 3
            sage: f(3/5)
            3/5
            sage: parent(f(3/5)) is K
            True
        """
        if S is QQ:
            return number_field_element_quadratic.Q_to_quadratic_field_element(self)
        else:
            return NumberField_absolute.coerce_map_from_impl(self, S)



    def discriminant(self, v=None):
        """
        Returns the discriminant of the ring of integers of the number field,
        or if v is specified, the determinant of the trace pairing
        on the elements of the list v.

        INPUT:
            v (optional) -- list of element of this number field
        OUTPUT:
            Integer if v is omitted, and Rational otherwise.

        EXAMPLES:
            sage: K.<i> = NumberField(x^2+1)
            sage: K.discriminant()
            -4
            sage: K.<a> = NumberField(x^2+5)
            sage: K.discriminant()
            -20
            sage: K.<a> = NumberField(x^2-5)
            sage: K.discriminant()
            5
        """
        if v is None:
            try:
                return self.__disc
            except AttributeError:
                d = self._D.squarefree_part()
                if d % 4 != 1:
                    d *= 4
                self.__disc = d
                return self.__disc
        else:
            return NumberField_generic.discriminant(self, v)

    def __reduce__(self):
        """
        This is used in pickling quadratic number fields.

        TESTS:
            sage: K.<z7> = QuadraticField(7)
            sage: L = loads(dumps(K))
            sage: L
            Number Field in z7 with defining polynomial x^2 - 7
            sage: L == K
            True
        """
        return NumberField_quadratic_v1, (self.polynomial(), self.variable_name())


    def is_galois(self):
        """
        Return True since all quadratic fields are automatically Galois.

        EXAMPLES:
            sage: QuadraticField(1234,'d').is_galois()
            True
        """
        return True

    def class_number(self, proof=None):
        r"""
        Return the size of the class group of self.

        If proof = False (*not* the default!) and the discriminant of the
        field is negative, then the following warning from the PARI
        manual applies: IMPORTANT WARNING: For $D<0$, this function
        may give incorrect results when the class group has a low
        exponent (has many cyclic factors), because implementing
        Shank's method in full generality slows it down immensely.

        EXAMPLES:
            sage: QuadraticField(-23,'a').class_number()
            3

        These are all the primes so that the class number of $\QQ(\sqrt{-p})$ is $1$:
            sage: [d for d in prime_range(2,300) if not is_square(d) and QuadraticField(-d,'a').class_number() == 1]
            [2, 3, 7, 11, 19, 43, 67, 163]

        It is an open problem to \emph{prove} that there are infinity
        many positive square-free $d$ such that $\QQ(\sqrt{d})$ has
        class number $1$:n
            sage: len([d for d in range(2,200) if not is_square(d) and QuadraticField(d,'a').class_number() == 1])
            121

        TESTS:
            sage: type(QuadraticField(-23,'a').class_number())
            <type 'sage.rings.integer.Integer'>
            sage: type(NumberField(x^3 + 23, 'a').class_number())
            <type 'sage.rings.integer.Integer'>
            sage: type(NumberField(x^3 + 23, 'a').extension(x^2 + 5, 'b').class_number())
            <type 'sage.rings.integer.Integer'>
            sage: type(CyclotomicField(10).class_number())
            <type 'sage.rings.integer.Integer'>
        """
        proof = proof_flag(proof)
        try:
            return self.__class_number
        except AttributeError:
            D = self.discriminant()
            if D < 0 and proof:
                self.__class_number = ZZ(pari("qfbclassno(%s,1)"%D))
            else:
                self.__class_number = ZZ(pari("qfbclassno(%s)"%D))
            return self.__class_number

    def hilbert_class_polynomial(self):
        r"""
        Returns a polynomial over $\QQ$ whose roots generate the
        Hilbert class field of this quadratic field.

        \note{Computed using PARI via Schertz's method.  This
        implementation is quite fast.}

        EXAMPLES:
            sage: K.<b> = QuadraticField(-23)
            sage: K.hilbert_class_polynomial()
            x^3 + x^2 - 1

            sage: K.<a> = QuadraticField(-431)
            sage: K.class_number()
            21
            sage: K.hilbert_class_polynomial()
            x^21 + x^20 - 13*x^19 - 50*x^18 + 592*x^17 - 2403*x^16 + 5969*x^15 - 10327*x^14 + 13253*x^13 - 12977*x^12 + 9066*x^11 - 2248*x^10 - 5523*x^9 + 11541*x^8 - 13570*x^7 + 11315*x^6 - 6750*x^5 + 2688*x^4 - 577*x^3 + 9*x^2 + 15*x + 1
        """
        f = pari('quadhilbert(%s))'%self.discriminant())
        g = QQ['x'](f)
        return g

    def hilbert_class_field(self, names):
        r"""
        Returns the Hilbert class field of this quadratic field as an
        absolute extension of $\QQ$.  For a polynomial that defines a
        relative extension see the \code{hilbert_class_polynomial}
        command.

        \note{Computed using PARI via Schertz's method.  This implementation
        is amazingly fast.}

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^2 + 23, 'a')
            sage: K.hilbert_class_polynomial()
            x^3 + x^2 - 1
            sage: K.hilbert_class_field('h')
            Number Field in h with defining polynomial x^6 + 2*x^5 + 70*x^4 + 90*x^3 + 1631*x^2 + 1196*x + 12743
        """
        f = self.hilbert_class_polynomial()
        C = self.composite_fields(NumberField(f,'x'),names)
        assert len(C) == 1
        return C[0]

def is_fundamental_discriminant(D):
    r"""
    Return True if the integer $D$ is a fundamental discriminant, i.e.,
    if $D \con 0,1\pmod{4}$, and $D\neq 0, 1$ and either (1) $D$ is square free
    or (2) we have $D\con 0\pmod{4}$ with $D/4 \con 2,3\pmod{4}$ and $D/4$
    square free.  These are exactly the discriminants of quadratic fields.

    EXAMPLES:
        sage: [D for D in range(-15,15) if is_fundamental_discriminant(D)]
        [-15, -11, -8, -7, -4, -3, 5, 8, 12, 13]
        sage: [D for D in range(-15,15) if not is_square(D) and QuadraticField(D,'a').disc() == D]
        [-15, -11, -8, -7, -4, -3, 5, 8, 12, 13]
    """
    d = D % 4
    if not (d in [0,1]):
        return False
    return D != 1 and  D != 0 and \
           (arith.is_squarefree(D) or \
            (d == 0 and (D//4)%4 in [2,3] and arith.is_squarefree(D//4)))


###################
# For pickling
###################


def NumberField_absolute_v1(poly, name, latex_name):
    """
    This is used in pickling generic number fields.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import NumberField_generic_v1
        sage: R.<x> = QQ[]
        sage: NumberField_generic_v1(x^2 + 1, 'i', 'i')
        Number Field in i with defining polynomial x^2 + 1
    """
    return NumberField_absolute(poly, name, latex_name, check=False)

NumberField_generic_v1 = NumberField_absolute_v1  # for historical reasons only (so old objects unpickle)

def NumberField_relative_v1(base_field, poly, name, latex_name):
    """
    This is used in pickling relative fields.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import NumberField_relative_v1
        sage: R.<x> = CyclotomicField(3)[]
        sage: NumberField_relative_v1(CyclotomicField(3), x^2 + 7, 'a', 'a')
        Number Field in a with defining polynomial x^2 + 7 over its base field
    """
    return NumberField_relative(base_field, poly, name, latex_name, check=False)

NumberField_extension_v1 = NumberField_relative_v1  # historical reasons only

def NumberField_cyclotomic_v1(zeta_order, name):
    """
    This is used in pickling cyclotomic fields.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import NumberField_cyclotomic_v1
        sage: NumberField_cyclotomic_v1(5,'a')
        Cyclotomic Field of order 5 and degree 4
        sage: NumberField_cyclotomic_v1(5,'a').variable_name()
        'a'
    """
    return NumberField_cyclotomic(zeta_order, name)

def NumberField_quadratic_v1(poly, name):
    """
    This is used in pickling quadratic fields.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import NumberField_quadratic_v1
        sage: R.<x> = QQ[]
        sage: NumberField_quadratic_v1(x^2 - 2, 'd')
        Number Field in d with defining polynomial x^2 - 2
    """
    return NumberField_quadratic(poly, name, check=False)


def put_natural_embedding_first(v):
    """
    Helper function for embeddings() functions for number fields.

    INPUT: a list of embeddings of a number field

    OUTPUT: None.  The list is altered in-place, so that, if possible,
            the first embedding has been switched with one of the
            others, so that if there is an embedding which preserves
            the generator names then it appears first.

    EXAMPLES:
        sage: K.<a> = CyclotomicField(7)
        sage: embs = K.embeddings(K)
        sage: [e(a) for e in embs] # already sorted
        [a, a^2, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1]
        sage: permuted_embs = [embs[i] for i in [1,2,3,4,5,0]]
        sage: [e(a) for e in permuted_embs] # natural map is not first
        [a^2, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1, a]
        sage: from sage.rings.number_field.number_field import put_natural_embedding_first
        sage: put_natural_embedding_first(permuted_embs)
        sage: [e(a) for e in permuted_embs] # now natural map is first
        [a, a^3, a^4, a^5, -a^5 - a^4 - a^3 - a^2 - a - 1, a^2]
    """
    for i in range(len(v)):
        phi = v[i]
        a = str(list(phi.domain().gens()))
        b = str(list(phi.im_gens()))
        if a == b:
            v[i] = v[0]
            v[0] = phi
            return



def refine_embedding(e, prec=None):
    """
    Given an embedding e: K->RR or CC, returns an equivalent embedding
    with higher precision.

    INPUT:
        e -- an embedding of a number field into either RR or CC (with
             some precision)
        prec -- (default None) the desired precision; if None, current
                precision is doubled.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field import refine_embedding
        sage: K = CyclotomicField(3)
        sage: e10 = K.complex_embedding(10)
        sage: e10.codomain().precision()
        10
        sage: e25 = refine_embedding(e10, prec=25)
        sage: e25.codomain().precision()
        25
    """
    K = e.domain()
    RC = e.codomain()
    prec_old = RC.precision()
    if prec is None:
        prec = 2*prec_old
    elif prec_old >= prec:
        return e

    # We first compute all the embeddings at the new precision:
    if sage.rings.real_mpfr.is_RealField(RC):
        elist = K.real_embeddings(prec)
    else:
        elist = K.complex_embeddings(prec)

    # Now we determine which is an extension of the old one; this
    # relies on the fact that coercing a high-precision root into a
    # field with lower precision will equal the lower-precision root!
    old_root = e(K.gen())
    diffs = [(RC(ee(K.gen()))-old_root).abs() for ee in elist]
    return elist[min(izip(diffs,count()))[1]]
