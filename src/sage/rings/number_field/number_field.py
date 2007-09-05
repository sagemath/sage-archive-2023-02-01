r"""
Number Fields

AUTHORS:
   -- William Stein (2004, 2005): initial version
   -- Steven Sivek (2006-05-12): added support for relative extensions
   -- William Stein (2007-09-04): major rewrite and documentation

This example follows one in the Magma reference manual:
    sage: K.<y> = NumberField(x^4 - 420*x^2 + 40000)
    sage: z = y^5/11; z
    420/11*y^3 - 40000/11*y
    sage: R.<y> = PolynomialRing(K)
    sage: f = y^2 + y + 1
    sage: L.<a> = K.extension(f); L
    Number Field in a with defining polynomial y^2 + y + 1 over its base field.
    sage: KL.<b> = NumberField([x^2 + x + 1, x^4 - 420*x^2 + 40000]); KL
    Number Field in b0 with defining polynomial x^4 + (-420)*x^2 + 40000 over its base field.
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

# There will be one running instance of GP for all
# number field calculations that use the interpreter.
from sage.interfaces.gp import Gp

import sage.libs.ntl.all as ntl
import sage.libs.pari.all as pari
import sage.interfaces.gap
import sage.misc.preparser
import sage.rings.arith
import sage.rings.complex_field
import sage.rings.ring
from sage.misc.latex import latex_variable_name, latex_varify

from class_group import ClassGroup
from galois_group import GaloisGroup

from sage.structure.element import is_Element

import sage.structure.parent_gens

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

from sage.structure.parent_gens import ParentWithGens
import number_field_element
from number_field_ideal import convert_from_zk_basis

from sage.libs.all import pari, pari_gen

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()

_nf_cache = {}
def NumberField(polynomial, name=None, check=True, names=None, all=False):
    r"""
    Return {\em the} number field defined by the given irreducible
    polynomial and with variable with the given name.  If check is
    True (the default), also verify that the defining polynomial is
    irreducible and over Q.

    INPUT:
        polynomial -- a polynomial over QQ or a number field, or
                      a list of polynomials.
        name -- a string (default: 'a'), the name of the generator
        check -- bool (default: True); do type checking and
                 irreducibility checking.
        all -- bool (default: False); if True, return a list
               of number fields, one for each factor.

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

    EXAMPLES: Constructing a relative number field
        sage: K.<a> = NumberField(x^2 - 2)
        sage: R.<t> = K[]
        sage: L = K.extension(t^3+t+a, 'b'); L
        Number Field in b with defining polynomial t^3 + t + a over its base field.
        sage: L.absolute_field()
        Number Field in b with defining polynomial x^6 + 2*x^4 + x^2 - 2
        sage: b = L.gen()
        sage: a*b
        -b^4 - b^2
        sage: L.lift_to_base(-3*b^3 - 3*b + 1)
        3*a + 1

    Constructing another number field:
        sage: k.<i> = NumberField(x^2 + 1)
        sage: R.<z> = k[]
        sage: m.<j> = NumberField(z^3 + i*z + 3)
        sage: m
        Number Field in j with defining polynomial z^3 + i*z + 3 over its base field.

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
        sage: M.<z> = NumberField([theta^2 + 3, theta^3 + 4]); M
        Number Field in z0 with defining polynomial theta^3 + 4 over its base field.
    """
    if name is None and names is None:
        raise TypeError, "You must specify the name of the generator."
    if not names is None:
        name = names

    name = sage.structure.parent_gens.normalize_names(1, name)

    if isinstance(polynomial, (list, tuple)):
        return MultiNumberField(polynomial, name, all=all)


    if not isinstance(polynomial, polynomial_element.Polynomial):
        try:
            polynomial = polynomial.polynomial(QQ)
        except (AttributeError, TypeError):
            raise TypeError, "polynomial (=%s) must be a polynomial."%repr(polynomial)

    if all:
        return [NumberField(f, name=name, check=check, names=names,all=False) for f, _ in polynomial.factor()]

    key = (polynomial, name)
    if _nf_cache.has_key(key):
        K = _nf_cache[key]()
        if not K is None: return K

    R = polynomial.base_ring()
    if R == ZZ:
        polynomial = QQ['x'](polynomial)
    elif isinstance(R, NumberField_generic):
        S = R.extension(polynomial, name)
        _nf_cache[key] = weakref.ref(S)
        return S

    if polynomial.degree() == 2:
        K = NumberField_quadratic(polynomial, name, check)
    else:
        K = NumberField_generic(polynomial, name, None, check)

    _nf_cache[key] = weakref.ref(K)
    return K


def MultiNumberField(v, names, all=False):
    """
    Return the number field defined by the polynomials or number
    fields in the list v.

    This is the field constructed first from v[0], then over that
    field from v[1], etc.  If all is False, then each v[i] must be
    irreducible over the previous fields.  Otherwise a list of all
    possible fields defined by all those polynomials is output.

    If names defines a variable name a, say, then the generators of
    the intermediate number fields are a0, a1, a2, ...

    INPUT:
        v -- a list of polynomials or number fields
        names -- variable name
        all -- bool (default: False), if True then return all possible
               number fields constructed from the given input data;
               otherwise just contstruct 1.

    OUTPUT:
        a single number field or a list of number fields

    EXAMPLES:
        sage: k.<a> = NumberField([x^2 + 1, x^2 + 3, x^2 + 5]); k
        Number Field in a1 with defining polynomial x^2 + 5 over its base field.

    Note -- because SAGE currently doesn't support relative extensions
    of relative extensions, the base field is an absolute field.
        sage: k.base_field()
        Number Field in a0 with defining polynomial x^4 + 8*x^2 + 4

    In the following examle the second polynomial reducible over the first, so
    we have to create all extensions or we get an error:
        sage: v = NumberField([x^3 - 2, x^3 - 2], names='a')
        Traceback (most recent call last):
        ...
        ValueError: defining polynomial (x^3 + -2) must be irreducible
        sage: v = NumberField([x^3 - 2, x^3 - 2], all=True, names='a'); v
        [Number Field in a1 with defining polynomial x + -a0 over its base field.,
         Number Field in a2 with defining polynomial x^2 + a0*x + a0^2 over its base field.]
        sage: v[0].absolute_field()
        Number Field in a1 with defining polynomial x^3 - 2
        sage: v[1].absolute_field()
        Number Field in a2 with defining polynomial x^6 + 108
        sage: v[1].absolute_field().galois_group()
        Galois group PARI group [6, -1, 2, "D_6(6) = [3]2"] of degree 6 of the number field Number Field in a2 with defining polynomial x^6 + 108


    We mix polynomial parent rings:
        sage: k.<y> = QQ[]
        sage: m = NumberField([y^3 + 2, x^2 + x + 1, y^3 - 3], 'beta')
        sage: m
        Number Field in beta1 with defining polynomial y^3 + -3 over its base field.
        sage: m.base_field ()
        Number Field in beta0 with defining polynomial y^6 + 3*y^5 + 6*y^4 + 3*y^3 + 9*y + 9

    """
    name = sage.structure.parent_gens.normalize_names(1, names)[0]
    if not isinstance(v, (list, tuple)):
        raise TypeError, "v must be a list or tuple"
    if len(v) == 0:
        return QQ
    if len(v) == 1:
        if all:
            f = v[0]
            if not isinstance(f, polynomial_element.Polynomial):
                f = QQ['x'](v[0])
            F = f.factor()
            return [NumberField(F[i][0], names=name+str(i)) for i in range(len(F))]
        else:
            return NumberField(v[0], names=names)
    f = v[-1]
    w = MultiNumberField(v[:-1], names=names, all=all)
    if is_NumberFieldExtension(w):
        w = w.absolute_field()
    if isinstance(f, polynomial_element.Polynomial):
        var = f.name()
    else:
        var = 'x'
    if all:
        R = w[-1][var]  # polynomial ring
        s = w[-1].variable_name()[len(name):]
    else:
        R = w[var]  # polynomial ring
        s = w.variable_name()[len(name):]

    f = R(f)
    i = 0
    if s:
        i = int(s) + 1
    else:
        i = 0

    if all:
        z = []
        for k in w:
            for g, _ in f.factor():
                z.append(k.extension(g, names=name+str(i)))
                i += 1
        return z
    else:
        return w.extension(f, name+str(i))


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
        sage: type(K)
        <class 'sage.rings.number_field.number_field.NumberField_quadratic'>
        sage: is_NumberField(K)
        True
    """
    D = QQ(D)
    if check:
        if D.is_square():
            raise ValueError, "D must not be a perfect square."
    R = polynomial_ring.PolynomialRing(QQ, 'x')
    f = R([-D, 0, 1])
    return NumberField(f, names, check=False)

def is_QuadraticField(x):
    r"""
    Return True if x is of the quadratic {\em number} field type.

    EXAMPLES:
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

def is_NumberFieldExtension(x):
    """
    Return True if x is an extension of a number field, i.e., a relative
    number field.

    EXAMPLES:
        sage: is_NumberFieldExtension(NumberField(x^2+1,'a'))
        False
        sage: k.<a> = NumberField(x^3 - 2)
        sage: l.<b> = k.extension(x^3 - 3); l
        Number Field in b with defining polynomial x^3 + -3 over its base field.
        sage: is_NumberFieldExtension(l)
        True
        sage: is_NumberFieldExtension(QQ)
        False
    """
    return isinstance(x, NumberField_extension)

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
        TypeError: no coercion of this rational to integer

    The degree must be positive.
        sage: CyclotomicField(0)
        Traceback (most recent call last):
        ...
        ValueError: n (=0) must be a positive integer

    The special case $n=1$ does \emph{not} return the rational numbers:
        sage: CyclotomicField(1)
        Cyclotomic Field of order 1 and degree 1
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

        If you use check=False, you avoid checking irreducibility of
        the defining polynomial, which can save time.
            sage: K.<a> = NumberField(x^2 - 1, check=False)

        It can also be dangerous:
            sage: (a-1)*(a+1)
            0
        """
        ParentWithGens.__init__(self, QQ, name)
        if not isinstance(polynomial, polynomial_element.Polynomial):
            raise TypeError, "polynomial (=%s) must be a polynomial"%repr(polynomial)

        if check:
            if not polynomial.is_irreducible():
                raise ValueError, "defining polynomial (%s) must be irreducible"%polynomial
            if not polynomial.parent().base_ring() == QQ:
                raise TypeError, "polynomial must be defined over rational field"
        if not polynomial.is_monic():
            raise NotImplementedError, "number fields for non-monic polynomials not yet implemented."

        self._assign_names(name)
        if latex_name is None:
            self.__latex_variable_name = latex_variable_name(self.variable_name())
        else:
            self.__latex_variable_name = latex_name
        self.__polynomial = polynomial
        self.__pari_bnf_certified = False
        self.__absolute_field = self

    def __reduce__(self):
        """
        TESTS:
            sage: Z = var('Z')
            sage: K.<w> = NumberField(Z^3 + Z + 1)
            sage: L = loads(dumps(K))
            sage: print L
            Number Field in w with defining polynomial Z^3 + Z + 1
            sage: print L == K
            True
        """
        return NumberField_generic_v1, (self.__polynomial, self.variable_name(), self.__latex_variable_name)

    def complex_embeddings(self, prec=53):
        r"""
        Return all homomorphisms of this ring into the approximate
        complex field with precision prec.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: f = x^5 + x + 17
            sage: k.<a> = NumberField(f)
            sage: v = k.complex_embeddings()
            sage: [phi(k.0^2) for phi in v]
            [0.921039066973047 - 3.07553311884578*I, 0.921039066973047 + 3.07553311884578*I, 2.97572074037668, -2.40889943716139 - 1.90254105303505*I, -2.40889943716139 + 1.90254105303505*I]
        """
        try:
            return self.__complex_embeddings[prec]
        except AttributeError:
            self.__complex_embeddings = {}
        except KeyError:
            pass
        CC = sage.rings.complex_field.ComplexField(prec)
        f = self.defining_polynomial().base_extend(CC)
        v = f.roots()
        e = [self.hom([a], check=False) for a in v]
        self.__complex_embeddings[prec] = e
        return e

    def latex_variable_name(self, name=None):
        """
        Return the latex representation of the variable name for
        this number field.

        EXAMPLES:

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
        """
        Return latex representation of this number field.  This is viewed
        as a polynomial quotient ring over a field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^13 - (2/3)*x + 3)
            sage: k._latex_()
            '\\mathbf{Q}[a]/(a^{13} - \\frac{2}{3}a + 3)'
            sage: latex(k)
            \mathbf{Q}[a]/(a^{13} - \frac{2}{3}a + 3)

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
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return number_field_element.NumberFieldElement(self, x.polynomial())
            return self._coerce_from_other_number_field(x)
        elif isinstance(x,str):
            return self._coerce_from_str(x)
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
            return number_field_element.NumberFieldElement(self, f[0])
        # todo: more general coercion if embedding have been asserted
        raise TypeError, "Cannot coerce %s into %s"%(x,self)

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

        Not any polynomial coerces in, e.g., not this one in characteristic 7.
            sage: f = GF(7)['y']([1,2,3]); f
            3*y^2 + 2*y + 1
            sage: K._coerce_non_number_field_element_in(f)
            Traceback (most recent call last):
            ...
            TypeError
        """
        if isinstance(x, (int, long, rational.Rational,
                              integer.Integer, pari_gen,
                              list)):
            return number_field_element.NumberFieldElement(self, x)

        try:
            if isinstance(x, polynomial_element.Polynomial):
                return number_field_element.NumberFieldElement(self, x)

            return number_field_element.NumberFieldElement(self, x._rational_())
        except (TypeError, AttributeError):
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
            (-7/8, <type 'sage.rings.number_field.number_field_element.NumberFieldElement'>)
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
            return number_field_element.NumberFieldElement(self, x)
        elif isinstance(x, number_field_element.NumberFieldElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return number_field_element.NumberFieldElement(self, x.list())
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
        """
        if not isinstance(other, NumberField_generic):
            return cmp(type(self), type(other))
        c = cmp(self.variable_name(), other.variable_name())
        if c: return c
        return cmp(self.__polynomial, other.__polynomial)

    def _ideal_class_(self):
        """
        Return the Python class used in defining ideals of the ring of
        integes of this number field.

        This function is required by the general ring/ideal machinery.

        EXAMPLES:
            sage: NumberField(x^2 + 2, 'c')._ideal_class_()
            <class 'sage.rings.number_field.number_field_ideal.NumberFieldIdeal'>
        """
        return sage.rings.number_field.number_field_ideal.NumberFieldIdeal

    def ideal(self, gens):
        r"""
        Return the ideal in $\mathcal{O}_K$ generated by gens.  This
        overrides the \code{sage.rings.ring.Field} method to use the
        \code{sage.rings.ring.Ring} one instead, since we're not really
        concerned with ideals in a field but in its ring of integers.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-2)
            sage: K.ideal([a])
            Fractional ideal (a) of Number Field in a with defining polynomial x^3 - 2
        """
        return sage.rings.ring.Ring.ideal(self, gens)

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

    def pari_polynomial(self):
        """
        PARI polynomial corresponding to polynomial that defines
        this field.   This is always a polynomial in the variable "x".

        EXAMPLES:
            sage: y = polygen(QQ)
            sage: k.<a> = NumberField(y^2 - 3/2*y + 5/3)
            sage: k.pari_polynomial()
            x^2 - 3/2*x + 5/3
        """
        try:
            return self.__pari_polynomial
        except AttributeError:
            self.__pari_polynomial = self.polynomial()._pari_()
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
            [x^2 + x + 1, [0, 1], -3, 1, [Mat([1, -0.50000000000000000000000000000000000000 + 0.86602540378443864676372317075293618347*I]), [1, 0.36602540378443864676372317075293618347; 1, -1.3660254037844386467637231707529361835], 0, [2, -1; -1, -1], [3, 2; 0, 1], [1, -1; -1, -2], [3, [2, -1; 1, 1]]], [-0.50000000000000000000000000000000000000 + 0.86602540378443864676372317075293618347*I], [1, x], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
            sage: pari(k)
            [x^2 + x + 1, [0, 1], -3, 1, [Mat([1, -0.50000000000000000000000000000000000000 + 0.86602540378443864676372317075293618347*I]), [1, 0.36602540378443864676372317075293618347; 1, -1.3660254037844386467637231707529361835], 0, [2, -1; -1, -1], [3, 2; 0, 1], [1, -1; -1, -2], [3, [2, -1; 1, 1]]], [-0.50000000000000000000000000000000000000 + 0.86602540378443864676372317075293618347*I], [1, x], [1, 0; 0, 1], [1, 0, 0, -1; 0, 1, 1, -1]]
        """
        if self.defining_polynomial().denominator() != 1:
            raise TypeError, "Unable to coerce number field defined by non-integral polynomial to PARI."
        return 'nfinit(%s)'%self.pari_polynomial()

    def pari_bnf(self, certify=False):
        """
        PARI big number field corresponding to this field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: len(k.pari_bnf())
            10
            sage: k.pari_bnf()[:4]
            [[;], matrix(0,7), [;], Mat([0.E-693 + 10.352073178770991947816442612760937457*I, 0.E-693 + 8.2487727536742446128967079885883377973*I, 0.E-693 + 9.9147352870231080237320951122610602744*I, 0.E-693 + 10.995574287564276334619251841478260095*I, 0.E-693 + 5.3558900891779742444967743036365769643*I, 0.E-693 + 4.3175978606849283409538655445296737394*I, 0.E-693 + 2.6516353273360649301184784208569512624*I])]
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
            if self.pari_bnf(certify=False).bnfcertify() != 1:
                raise ValueError, "The result is not correct according to bnfcertify"
            self.__pari_bnf_certified = True
        return self.__pari_bnf_certified

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

    def class_group(self, proof=True, names='c'):
        r"""
        Return the class group of this field.

        INPUT:
            proof -- if True (which is the default), then compute
                       the classgroup provably correctly.
            names -- names of the generators of this class group.

        OUTPUT:
            The class group of this number field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 23); k
            Number Field in a with defining polynomial x^2 + 23
            sage: G = k.class_group(); G
            Multiplicative Abelian Group isomorphic to C3 as the class group of Number Field in a with defining polynomial x^2 + 23
            sage: G.number_field()
            Number Field in a with defining polynomial x^2 + 23
            sage: G is k.class_group()
            True
            sage: G is k.class_group(proof=False)
            False
            sage: G.gens()
            (c,)

        There can be multiple generators:
            sage: k.<a> = NumberField(x^2 + 20072)
            sage: G = k.class_group(); G
            Multiplicative Abelian Group isomorphic to C2 x C2 x C19 as the class group of Number Field in a with defining polynomial x^2 + 20072
            sage: G.gens()
            (c0, c1)

        You can name the generators during construction:
            sage: G.<Z0,Z1> = k.class_group(); G.gens()
            (Z0, Z1)

        Assigning generator names without having to know how many
        there will be:
            sage: k.class_group(names='W').gens()
            (W0, W1)

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
        try:
            return self.__class_group[proof, names]
        except KeyError:
            pass
        except AttributeError:
            self.__class_group = {}
        k = self.pari_bnf(proof)
        s = str(k.getattr('clgp'))
        s = s.replace(";",",")
        s = eval(s)
        G = ClassGroup(s[1], names, self)
        self.__class_group[proof,names] = G
        return G

    def class_number(self, proof=True):
        """
        Return the class number of this number field, as an integer.

        INPUT:
            proof -- bool (default: True)

        EXAMPLES:
            sage: NumberField(x^2 + 23, 'a').class_number()
            3
            sage: NumberField(x^2 + 163, 'a').class_number()
            1
            sage: NumberField(x^3 + x^2 + 997*x + 1, 'a').class_number(proof=False)
            1539
        """
        return self.class_group(proof).order()

    def composite_fields(self, other, names):
        """
        List of all possible composite number fields formed from self
        and other.

        INPUT:
            other -- a number field
            names -- generator name for composite fields

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 + 2)
            sage: m.<b> = NumberField(x^3 + 2)
            sage: k.composite_fields(m, 'c')
            [Number Field in c with defining polynomial x^3 - 2,
             Number Field in c with defining polynomial x^6 - 40*x^3 + 1372]
        """
        if not isinstance(other, NumberField_generic):
            raise TypeError, "other must be a number field."
        f = self.pari_polynomial()
        g = other.pari_polynomial()
        C = f.polcompositum(g)
        R = self.polynomial().parent()
        C = [R(h) for h in C]
        return [NumberField(h, names) for h in C]

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
            sage: d = k.different(); d
            Fractional ideal (a) of Number Field in a with defining polynomial x^2 + 23
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
            Fractional ideal (2*b) of Number Field in b with defining polynomial x^2 - 123
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
                self.__disc = QQ(str(self.pari_nf()[2]))
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

    def elements_of_norm(self, n, proof=True):
        r"""
        Return a list of solutions modulo units of positive norm to
        $Norm(a) = n$, where a can be any integer in this number field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2+1)
            sage: K.elements_of_norm(3)
            []
            sage: K.elements_of_norm(50)
            [7*a - 1, -5*a + 5, a - 7]           # 32-bit
            [7*a - 1, -5*a + 5, -7*a - 1]        # 64-bit
        """
        B = self.pari_bnf(proof).bnfisintnorm(n)
        R = self.polynomial().parent()
        return [self(QQ['x'](R(g))) for g in B]

    def extension(self, poly, name=None, names=None):
        """
        Return the relative extension of this field by a given polynomial.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 - 2)
            sage: R.<t> = K[]
            sage: L.<b> = K.extension(t^2 + a); L
            Number Field in b with defining polynomial t^2 + a over its base field.

        We create another extension.
            sage: k.<a> = NumberField(x^2 + 1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: y = var('y')
            sage: m.<b> = k.extension(y^2 + 2); m
            Number Field in b with defining polynomial y^2 + 2 over its base field.
            sage: b.minpoly()
            x^4 + 6*x^2 + 1
            sage: b.minpoly('z')
            z^4 + 6*z^2 + 1
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
        return NumberField_extension(self, poly, str(name))

    def factor_integer(self, n):
        r"""
        Ideal factorization of the principal ideal of the ring
        of integers generated by $n$.

	EXAMPLE:
        Here we show how to factor gaussian integers.
        First we form a number field defined by $x^2 + 1$:

            sage: K.<I> = NumberField(x^2 + 1); K
            Number Field in I with defining polynomial x^2 + 1

        Here are the factors:

	    sage: fi, fj = K.factor_integer(13); fi,fj
            ((Fractional ideal (3*I - 2) of Number Field in I with defining polynomial x^2 + 1, 1),
            (Fractional ideal (-3*I - 2) of Number Field in I with defining polynomial x^2 + 1, 1))

        Now we extract the reduced form of the generators:

	    sage: zi = fi[0].gens_reduced()[0]; zi
            3*I - 2
	    sage: zj = fj[0].gens_reduced()[0]; zj
            -3*I - 2

        We recover the integer that was factor in $\Z[i]$

	    sage: zi*zj
            13

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
            self.__gen = number_field_element.NumberFieldElement(self, X)
            return self.__gen

    def is_field(self):
        """
        Return True since a number field is a field.

        EXAMPLES:
            sage: NumberField(x^5 + x + 3, 'c').is_field()
            True
        """
        return True

    def galois_group(self, pari_group = True, use_kash=False):
        r"""
        Return the Galois group of the Galois closure of this number
        field as an abstract group.

        For more (important!) documentation, so the documentation
        for Galois groups of polynomials over $\Q$, e.g., by
        typing \code{K.polynomial().galois_group?}, where $K$
        is a number field.

        EXAMPLES:
            sage: k.<b> = NumberField(x^2 - 14)
            sage: k.galois_group ()
            Galois group PARI group [2, -1, 1, "S2"] of degree 2 of the number field Number Field in b with defining polynomial x^2 - 14

            sage: NumberField(x^3-2, 'a').galois_group(pari_group=True)
            Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the number field Number Field in a with defining polynomial x^3 - 2

            sage: NumberField(x-1, 'a').galois_group(pari_group=False)    # optional database_gap package
            Galois group Transitive group number 1 of degree 1 of the number field Number Field in a with defining polynomial x - 1
            sage: NumberField(x^2+2, 'a').galois_group(pari_group=False)  # optional database_gap package
            Galois group Transitive group number 1 of degree 2 of the number field Number Field in a with defining polynomial x^2 + 2
            sage: NumberField(x^3-2, 'a').galois_group(pari_group=False)  # optional database_gap package
            Galois group Transitive group number 2 of degree 3 of the number field Number Field in a with defining polynomial x^3 - 2
        """
        try:
            return self.__galois_group[pari_group, use_kash]
        except KeyError:
            pass
        except AttributeError:
            self.__galois_group = {}

        G = self.polynomial().galois_group(pari_group = pari_group, use_kash = use_kash)
        H = GaloisGroup(G, self)
        self.__galois_group[pari_group, use_kash] = H
        return H


    def integral_basis(self):
        """
        Return a list of elements of this number field that are a basis
        for the full ring of integers.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 + 10*x + 1)
            sage: K.integral_basis()
            [1, a, a^2, a^3, a^4]

        Next we compute the ring of integers of a cubic field in which 2
        is an "essential discriminant divisor", so the ring of integers
        is not generated by a single element.
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: K.integral_basis()
            [1, a, 1/2*a^2 + 1/2*a]
        """
        try:
            return self.__integral_basis
        except AttributeError:
            f = self.pari_polynomial()
            B = f.nfbasis()
            R = self.polynomial().parent()
            self.__integral_basis = [self(R(g).list()) for g in B]
        return self.__integral_basis

    def narrow_class_group(self, proof = True):
        r"""
        Return the narrow class group of this field.

        EXAMPLES:
            sage: NumberField(x^3+x+9, 'a').narrow_class_group()
            Multiplicative Abelian Group isomorphic to C2
        """
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
            sage: M.<a> = NumberField([y^2 +1 , y^3 + 97]); M
            Number Field in a0 with defining polynomial y^3 + 97 over its base field.
            sage: M.polynomial_ring()
            Univariate Polynomial Ring in y over Number Field in a with defining polynomial y^2 + 1
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

    def regulator(self, proof=True):
        """
        Return the regulator of this number field.

        Note that PARI computes the regulator to higher precision than
        the SAGE default.

        EXAMPLES:
            sage: NumberField(x^2-2, 'a').regulator()
            0.88137358701954305
            sage: NumberField(x^4+x^3+x^2+x+1, 'a').regulator()
            0.96242365011920694
        """
        try:
            return self.__regulator
        except AttributeError:
            k = self.pari_bnf(proof)
            s = str(k.getattr('reg'))
            self.__regulator = eval(s)
        return self.__regulator

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.

        EXAMPLES:
            sage: NumberField(x^2+1, 'a').signature()
            (0, 1)
            sage: NumberField(x^3-2, 'a').signature()
            (1, 1)
            sage: CyclotomicField(7).signature()
            (0, 3)
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

    def units(self, proof = True):
        """
        Return generators for the unit group modulo torsion.

        ALGORITHM: Uses PARI's bnfunit command.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: A = x^4 - 10*x^3 + 20*5*x^2 - 15*5^2*x + 11*5^3
            sage: K = NumberField(A, 'a')
            sage: K.units()
            [8/275*a^3 - 12/55*a^2 + 15/11*a - 2]
        """
        try:
            return self.__units
        except AttributeError:
            B = self.pari_bnf(proof).bnfunit()
            R = self.polynomial().parent()
            self.__units = [self(R(g)) for g in B]
            return self.__units


    def zeta(self, n=2, all=False):
        """
        Return an n-th root of unity in this field.  If all is True,
        return all of them.

        INPUT:
            n -- positive integer
            all -- bool, default: False.  If True, return a list
                   of all n-th roots of 1)

        If there are no n-th roots of unity in self (and all is
        False), this function raises an ArithmeticError exception.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^2 + 3, 'zeta3')
            sage: K.zeta(1)
            1
            sage: K.zeta(2)
            -1
            sage: K.zeta(2, all=True)
            [-1]
            sage: K.zeta(3)
            1/2*zeta3 - 1/2
            sage: K.zeta(3, all=True)
            [1/2*zeta3 - 1/2, -1/2*zeta3 - 1/2]
            sage: K.zeta(4)
            Traceback (most recent call last):
            ...
            ArithmeticError: There are no 4-th roots of unity self.

            sage: r.<x> = QQ[]
            sage: K.<a> = NumberField(x^2+1)
            sage: K.zeta(4)
            a
            sage: K.zeta(4,all=True)
            [a, -a]
            sage: K.zeta(3)
            Traceback (most recent call last):
            ...
            ArithmeticError: There are no 3-th roots of unity self.
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
            field = self.__absolute_field
            f = field.polynomial_ring().cyclotomic_polynomial(n)
            F = polynomial_ring.PolynomialRing(field, 'x')(f)
            R = F.roots()
            if len(R) == 0:
                if all:
                    return []
                else:
                    raise ArithmeticError, "There are no %s-th roots of unity self."%n
            if all:
                return [r[0] for r in R]
            else:
                return R[0][0]

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



class NumberField_extension(NumberField_generic):
    """
    EXAMPLES:
        sage: K.<a> = NumberField(x^3 - 2)
        sage: t = K['x'].gen()
        sage: L.<b> = K.extension(t^2+t+a); L
        Number Field in b with defining polynomial x^2 + x + a over its base field.
    """
    def __init__(self, base, polynomial, name,
                 latex_name=None, names=None, check=True):
        """
        Note: polynomial must be defined in the ring \code{K['x']}, where
        K is the base field.
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
            raise ValueError, "The polynomial must be defined over the base field"

        if check:
            if not polynomial.is_irreducible():
                raise ValueError, "defining polynomial (%s) must be irreducible"%polynomial

        # Generate the nf and bnf corresponding to the base field
        # defined as polynomials in y, e.g. for rnfisfree

        # Convert the polynomial defining the base field into a
        # polynomial in y to satisfy PARI's ordering requirements.
        # NOTE: This might not work properly if the base field is not
        #       defined by a polynomial in one variable.  But currently
        #       they are all defined in one variable, so no problem!

        Qx = base.polynomial().parent()
        Qy = (base.polynomial().base_ring())['y']
        phi = Qx.hom([Qy.gen()])
        base_polynomial_y = phi(base.polynomial())

        self.__base_nf = pari(base_polynomial_y).nfinit()
        self.__base_bnf = pari(base_polynomial_y).bnfinit()

        # Use similar methods to convert the polynomial defining the
        # relative extension into a polynomial in x, with y denoting
        # the generator of the base field.
        # NOTE: This should be rewritten if there is a way to extend
        #       homomorphisms K -> K' to homomorphisms K[x] -> K'[x].
        base_field_y = NumberField(base.polynomial(), 'y')
        Kx = base_field_y['x']
        i = base.hom([base_field_y.gen()]) # inclusion K -> K' with a -> y
        rel_coeffs = [i(c) for c in polynomial.coeffs()]
        polynomial_y = Kx(rel_coeffs)

        self.__pari_relative_polynomial = pari(str(polynomial_y))
        self.__rnf = self.__base_nf.rnfinit(self.__pari_relative_polynomial)

        self.__base_field = base
        NumberField_generic.__init__(self, self.absolute_polynomial(), name=name, latex_name=latex_name, check=False)

        self._assign_names(name)
        self.__relative_polynomial = polynomial
        self.__pari_bnf_certified = False

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
        return NumberField_extension_v1, (self.__base_field, self.polynomial(), self.variable_name(),
                                          self.latex_variable_name())

    def _repr_(self):
        return "Number Field in %s with defining polynomial %s over its base field."%(self.variable_name(), self.polynomial())

        #return "Extension by %s of the Number Field in %s with defining polynomial %s"%(
        #self.polynomial(), self.base_field().variable_name(),
        #    self.base_field().polynomial())

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
        Coerce x into this number field.
        """
        if isinstance(x, number_field_element.NumberFieldElement):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return number_field_element.NumberFieldElement(self, x.polynomial())
            if x.parent() == self.base_field():
                return self.__base_inclusion(x)

        if not isinstance(x, (int, long, rational.Rational,
                              integer.Integer, pari_gen,
                              polynomial_element.Polynomial,
                              list)):
            return self.base_field()(x)
            #raise TypeError, "Cannot coerce %s into %s"%(x,self)

        return number_field_element.NumberFieldElement(self, x)

    def _coerce_impl(self, x):
        if isinstance(x, number_field_element.NumberFieldElement):
            if x.parent() == self:
                return x
            if x.parent() == self.base_field():
                return self.__base_inclusion(x)
        elif isinstance(x, (rational.Rational, integer.Integer, int, long)):
            return number_field_element.NumberFieldElement(self, x)
        raise TypeError

    def __base_inclusion(self, element):
        """
        Given an element of the base field, give its inclusion into this
        extension (according to PARI's rnfeltreltoabs) in terms of the
        generator of this field.
        """
        if not number_field_element.is_NumberFieldElement(element):
            raise TypeError, "element must be a NumberFieldElement"
        if element.parent() != self.base_field():
            raise TypeError, "element must belong to the base field"
        base_field_y = NumberField(self.base_field().polynomial(), 'y')
        phi = self.base_field().hom([base_field_y.gen()])
        expr_x = self.pari_rnf().rnfeltreltoabs(str(phi(element)))

        # Convert to a polynomial in x, then to one in gen(), and return it
        return self(QQ['x'](str(expr_x).replace('^','**')))

    def _ideal_class_(self):
        return sage.rings.number_field.number_field_ideal.NumberFieldIdeal_rel

    def _pari_base_bnf(self, proof=False):
        # No need to certify the same field twice, so we'll just check
        # that the base field is certified.
        if proof:
            self.base_field().pari_bnf_certify()
        return self.__base_bnf

    def _pari_base_nf(self):
        return self.__base_nf

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "Only one generator."
        try:
            return self.__gen
        except AttributeError:
            X = rational_field.RationalField()['x'].gen()
            self.__gen = number_field_element.NumberFieldElement(self, X)
            return self.__gen

    def gen_relative(self):
        """
        Return root of defining polynomial, which is a generator of
        the relative number field over the base.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2+1); k
            Number Field in a with defining polynomial x^2 + 1
            sage: y = polygen(k)
            sage: m.<b> = k.extension(y^2+3); m
            Number Field in b with defining polynomial x^2 + 3 over its base field.
            sage: c = m.gen_relative(); c
            1/4*b^3 + 5/2*b
            sage: c^2 + 3
            0
            sage: m.gen()
            b
        """
        try:
            return self.__gen_relative
        except AttributeError:
            rnf = self.pari_rnf()
            f = (pari('x') - rnf[10][2]*rnf[10][1]).lift()
            self.__gen_relative = number_field_element.NumberFieldElement(self, f)
            return self.__gen_relative

    def pari_polynomial(self):
        """
        PARI polynomial corresponding to polynomial that defines
        this field.
        """
        try:
            return self.__pari_polynomial
        except AttributeError:
            self.__pari_polynomial = self.absolute_polynomial()._pari_()
            return self.__pari_polynomial

    def pari_rnf(self):
        return self.__rnf

    def pari_relative_polynomial(self):
        return self.__pari_relative_polynomial

    def absolute_field(self, name=None):
        r"""
        Return this field as an extension of $\Q$ rather than an
        extension of the base field.
        """
        try:
            return self.__absolute_field
        except AttributeError:
            if name is None:
                name = self.variable_name()
            self.__absolute_field = NumberField(self.absolute_polynomial(), name)
            return self.__absolute_field

    def absolute_polynomial(self):
        r"""
        Return the polynomial over $\Q$ which defines this field as an
        extension of the rational numbers.
        """
        try:
            return self.__absolute_polynomial
        except AttributeError:
            pbn = self._pari_base_nf()
            prp = self.pari_relative_polynomial()
            pari_poly = pbn.rnfequation(prp)
            R = self.base_field().polynomial().parent()
            self.__absolute_polynomial = R(pari_poly)
            return self.__absolute_polynomial

    def base_field(self):
        return self.__base_field

    def base_ring(self):
        return self.base_field()

    def discriminant(self, proof=True):
        """
        Return the relative discriminant of this extension $L/K$ as
        an ideal of $K$.  If you want the (rational) discriminant of
        $L/Q$, use e.g. \code{L.absolute_field().discriminant()}.

        Note that this uses PARI's \code{rnfdisc} function, which
        according to the documentation takes an \code{nf} parameter in
        GP but a \code{bnf} parameter in the C library.  If the C
        library actually accepts an \code{nf}, then this function
        should be fixed and the \code{proof} parameter removed.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: K.<i> = NumberField(x^2+1)
            sage: t = K['x'].gen()
            sage: L.<b> = K.extension(t^4-i)
            sage: L.discriminant()
            Fractional ideal (256) of Number Field in i with defining polynomial x^2 + 1
        """
        bnf = self._pari_base_bnf(proof)
        K = self.base_field()
        R = K.polynomial().parent()
        D, d = bnf.rnfdisc(self.pari_relative_polynomial())
        return K.ideal([ K(R(x)) for x in convert_from_zk_basis(K, D) ])

    disc = discriminant

    def extension(self, poly, name='b'):
        """
        Raise a NotImplemented error, since relative extensions of relative
        extensions are not yet supported.
        """
        raise NotImplementedError, "relative extensions of relative extensions are not supported"

    def galois_group(self, pari_group = True, use_kash=False):
        r"""
        Return the Galois group of the Galois closure of this number
        field as an abstract group.  Note that even though this is an
        extension $L/K$, the group will be computed as if it were $L/\Q$.

        For more (important!) documentation, so the documentation
        for Galois groups of polynomials over $\Q$, e.g., by
        typing \code{K.polynomial().galois_group?}, where $K$
        is a number field.

        EXAMPLE:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^2 + 1)
            sage: R.<t> = PolynomialRing(K)
            sage: L = K.extension(t^5-t+a, 'b')
            sage: L.galois_group()
            Galois group PARI group [240, -1, 22, "S(5)[x]2"] of degree 10 of the number field Number Field in b with defining polynomial t^5 + (-1)*t + a over its base field.
        """
        try:
            return self.__galois_group[pari_group, use_kash]
        except KeyError:
            pass
        except AttributeError:
            self.__galois_group = {}

        G = self.absolute_polynomial().galois_group(pari_group = pari_group,
                                                    use_kash = use_kash)
        H = GaloisGroup(G, self)
        self.__galois_group[pari_group, use_kash] = H
        return H


    def is_free(self, proof=True):
        r"""
        Determine whether or not $L/K$ is free (i.e. if $\mathcal{O}_L$ is
        a free $\mathcal{O}_K$-module).

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^2+6)
            sage: L.<b> = K.extension(K['x'].gen()^2 + 3)    ## extend by x^2+3
            sage: L.is_free()
            False
        """
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
        if str(poly_xy).find('x') >= 0:
            raise ValueError, "The element %s is not in the base field"%element
        return self.base_field()( QQ['y'](poly_xy) )

    def polynomial(self):
        return self.__relative_polynomial



class NumberField_cyclotomic(NumberField_generic):
    """
    Create a cyclotomic extension of the rational field.

    The command CyclotomicField(n) creates the n-th cyclotomic
    field, got by adjoing an n-th root of unity to the rational
    field.

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
        sage: print z6*x^3/(z6 + x)
        zeta12^2*x^3/(x + zeta12^2)
    """
    def __init__(self, n, names):
        f = QQ['x'].cyclotomic_polynomial(n)
        if names[0][:4] == 'zeta':
            latex_name = "\\zeta_{%s}"%n
        else:
            latex_name = None
        NumberField_generic.__init__(self, f,
                                     name= names,
                                     latex_name=latex_name,
                                     check=False)
        n = integer.Integer(n)
        zeta = self.gen()
        zeta._set_multiplicative_order(n)
        self.__zeta_order = n

    def __reduce__(self):
        """
        TESTS:
            sage: K.<zeta7> = CyclotomicField(7)
            sage: L = loads(dumps(K))
            sage: print L
            Cyclotomic Field of order 7 and degree 6
            sage: print L == K
            True
        """
        return NumberField_cyclotomic_v1, (self.__zeta_order, self.variable_name())

    def _repr_(self):
        return "Cyclotomic Field of order %s and degree %s"%(
                self.zeta_order(), self.degree())

    def _latex_(self):
        return "%s(\\zeta_{%s})"%(latex(QQ), self.__zeta_order)

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

        Coercion of GAP cyclotomic elements is also fully supported.


        """
        if isinstance(x, number_field_element.NumberFieldElement):
            if isinstance(x.parent(), NumberField_cyclotomic):
                return self._coerce_from_other_cyclotomic_field(x)
            else:
                return self._coerce_from_other_number_field(x)
        elif sage.interfaces.gap.is_GapElement(x):
            return self._coerce_from_gap(x)
        elif isinstance(x,str):
            return self._coerce_from_str(x)
        else:
            return self._coerce_non_number_field_element_in(x)

    def _coerce_from_other_cyclotomic_field(self, x, only_canonical=False):
        """
        Coerce an element x of a cyclotomic field into self, if at all possible.

        INPUT:
            x -- number field element
            only_canonical -- bool (default: False); Attempt to work, even in some
                   cases when x is not in a subfield of the cyclotomics (as long as x is
                   a root of unity).
        """
        K = x.parent()
        if K is self:
            return x
        elif K == self:
            return number_field_element.NumberFieldElement(self, x.polynomial())
        n = K.zeta_order()
        m = self.zeta_order()
        if m % n == 0:   # easy case
            # pass this off to a method in the element class
            # it can be done very quickly and easily by the pyrex<->NTL interface there
            return x._lift_cyclotomic_element(self)
        else:
            if only_canonical:
                raise TypeError
            n = x.multiplicative_order()
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
                y = K(self.zeta())
                z = y
                for r in xrange(y.multiplicative_order()):
                    if z == x:
                        return self.zeta()**(r+1)
                    z *= y
            raise TypeError, "Cannot coerce %s into %s"%(x,self)
        return number_field_element.NumberFieldElement(self, g)

    def _coerce_from_gap(self, x):
        """
        Attempt to coerce a GAP number field element into this cyclotomic field.
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
        Canonical coercion of x into self.

        Elements of other compatible cyclotomic fields coerce in, as do elements
        of the rings that coerce to all number fields (e.g., integers, rationals).
        """
        if isinstance(x, number_field_element.NumberFieldElement) and \
                isinstance(x.parent(), NumberField_cyclotomic):
            return self._coerce_from_other_cyclotomic_field(x, only_canonical=True)
        return NumberField_generic._coerce_impl(self, x)

    def complex_embedding(self, prec=53):
        r"""
        Return the embedding of this cyclotomic field into the
        approximate complex field with precision prec obtained by
        sending the generator $\zeta$ of self to exp(2*pi*i/n), where
        $n$ is the multiplicative order of $\zeta$.

        EXAMPLES:
            sage: C = CyclotomicField(4)
            sage: C.complex_embedding()
            Ring morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Complex Field with 53 bits of precision
              Defn: zeta4 |--> 6.12323399573677e-17 + 1.00000000000000*I

        Note in the example above that the way zeta is computed (using
        sin and cosine in MPFR) means that only the prec bits of the
        number after the decimal point are valid.

            sage: K = CyclotomicField(3)
            sage: phi = K.complex_embedding (10)
            sage: phi(K.0)
            -0.50 + 0.87*I
            sage: phi(K.0^3)
            1.0
            sage: phi(K.0^3 - 1)
            0
            sage: phi(K.0^3 + 7)
            8.0
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        return self.hom([CC.zeta(self.zeta_order())], check=False)

    def complex_embeddings(self, prec=53):
        r"""
        Return all embeddings of this cyclotomic field into the
        approximate complex field with precision prec.

        EXAMPLES:
            sage: C = CyclotomicField(4)
            sage: C.complex_embeddings()
            [Ring morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Complex Field with 53 bits of precision
              Defn: zeta4 |--> 6.12323399573677e-17 + 1.00000000000000*I, Ring morphism:
              From: Cyclotomic Field of order 4 and degree 2
              To:   Complex Field with 53 bits of precision
              Defn: zeta4 |--> -0.000000000000000183697019872103 - 1.00000000000000*I]
        """
        CC = sage.rings.complex_field.ComplexField(prec)
        n = self.zeta_order()
        z = CC.zeta(self.zeta_order())
        X = [m for m in range(n) if sage.rings.arith.gcd(m,n) == 1]
        return [self.hom([z**n], check=False) for n in X]

    def next_split_prime(self, p=2):
        """
        Return the next prime integer $p$ that splits completely in
        this cyclotomic field (and does not ramify).

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: K.next_split_prime(7)
            13
        """
        n = self.zeta_order()
        while True:
            p = sage.rings.arith.next_prime(p)
            if p % n == 1:
                return p

    def integral_basis(self):
        """
        Return a list of elements of this number field that are a basis
        for the full ring of integers.
        """
        try:
            return self.__integral_basis
        except AttributeError:
            z = self.gen()
            a = self(1)
            B = []
            for n in xrange(self.degree()):
                B.append(a)
                a *= z
            self.__integral_basis = B
        return self.__integral_basis


    def zeta_order(self):
        return self.__zeta_order

    def multiplicative_order_table(self):
        try:
            return self.__multiplicative_order_table
        except AttributeError:
            t = {}
            x = self(1)
            n = self.zeta_order()
            m = 0
            zeta = self.zeta()
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
        number field, if there is one.  Raises a ValueError if there
        is not.

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

class NumberField_quadratic(NumberField_generic):
    """
    Create a quadratic extension of the rational field.

    The command QuadraticExtension(a) creates the field Q(sqrt(a)).

    EXAMPLES:
        sage: QuadraticField(3, 'a')
        Number Field in a with defining polynomial x^2 - 3
        sage: QuadraticField(-4, 'b')
        Number Field in b with defining polynomial x^2 + 4
    """
    def __init__(self, polynomial, name=None, check=True):
        NumberField_generic.__init__(self, polynomial, name=name, check=check)

    def __reduce__(self):
        """
        TESTS:
            sage: K.<z7> = QuadraticField(7)
            sage: L = loads(dumps(K))
            sage: print L
            Number Field in z7 with defining polynomial x^2 - 7
            sage: print L == K
            True
        """
        return NumberField_quadratic_v1, (self.polynomial(), self.variable_name())

    def class_number(self, proof = True):
        """
        Return the size of the class group of self.

        If proof = False (not the default) and the discriminant of the
        field is negative, then the following warning from the PARI
        manual applies: IMPORTANT WARNING: For D<0, this function may
        give incorrect results when the class group has a low exponent
        (has many cyclic factors), because implementing Shank's method
        in full generality slows it down immensely.
        """
        try:
            return self.__class_number
        except AttributeError:
            D = self.discriminant()
            if D < 0 and proof:
                self.__class_number = pari("qfbclassno(%s,1)"%D).python()
            else:
                self.__class_number = pari("qfbclassno(%s)"%D).python()
            return self.__class_number

    def hilbert_class_polynomial(self):
        r"""
        Returns a polynomial over $\Q$ whose roots generate the
        Hilbert class field of this quadratic field.

        \note{Computed using PARI via Schertz's method.  This
        implementation is quite fast.}

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: K = NumberField(x^2 + 23, 'a')
            sage: K.hilbert_class_polynomial()
            x^3 + x^2 - 1

            sage: K = NumberField(x^2 + 431, 'a')
            sage: K.hilbert_class_polynomial()
            x^21 + x^20 - 13*x^19 - 50*x^18 + 592*x^17 - 2403*x^16 + 5969*x^15 - 10327*x^14 + 13253*x^13 - 12977*x^12 + 9066*x^11 - 2248*x^10 - 5523*x^9 + 11541*x^8 - 13570*x^7 + 11315*x^6 - 6750*x^5 + 2688*x^4 - 577*x^3 + 9*x^2 + 15*x + 1
        """
        f = pari('quadhilbert(%s))'%self.discriminant())
        g = QQ['x'](f)
        return g

    def hilbert_class_field(self, names):
        r"""
        Returns the Hilbert class field of this quadratic
        field as an absolute extension of $\Q$.  For a polynomial
        that defines a relative extension see the
        \code{hilbert_class_polynomial} command.

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
    d = D % 4
    if not (d in [0,1]):
        return False
    return D != 1 and  D != 0 and \
           (arith.is_squarefree(D) or \
            (d == 0 and (D//4)%4 in [2,3] and arith.is_squarefree(D//4)))


###################
# For pickling
###################

def NumberField_generic_v1(poly, name, latex_name):
    return NumberField_generic(poly, name, latex_name, check=False)

def NumberField_extension_v1(base_field, poly, name, latex_name):
    return NumberField_extension(base_field, poly, name, latex_name, check=False)

def NumberField_cyclotomic_v1(zeta_order, name):
    return NumberField_cyclotomic(zeta_order, name)

def NumberField_quadratic_v1(poly, name):
    return NumberField_quadratic(poly, name, check=False)
