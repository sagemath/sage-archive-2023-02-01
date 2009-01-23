r"""
Relative Number Fields

AUTHORS:
   -- William Stein (2004, 2005): initial version
   -- Steven Sivek (2006-05-12): added support for relative extensions
   -- William Stein (2007-09-04): major rewrite and documentation
   -- Robert Bradshaw (2008-10): specified embeddings into ambient fields
   -- Nick Alexander (2009-01): modernize coercion implementation

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

#*****************************************************************************
#       Copyright (C) 2004-2009 William Stein <wstein@gmail.com>
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
from sage.categories.map import Map

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
import sage.rings.real_lazy

from sage.rings.integer_mod import mod

import sage.rings.ring
from sage.misc.latex import latex_variable_name, latex_varify

from class_group import ClassGroup
from galois_group import GaloisGroup

from sage.structure.element import is_Element
from sage.categories.map import is_Map
from sage.structure.sequence import Sequence

import sage.structure.parent_gens

from sage.structure.proof.proof import get_flag
import maps
import number_field_morphisms
from itertools import count, izip

from sage.rings.integer_ring import IntegerRing

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
from sage.rings.number_field.number_field import NumberField, NumberField_generic, put_natural_embedding_first, proof_flag
from sage.rings.number_field.number_field_base import is_NumberField

from sage.rings.number_field.number_field_ideal_rel import NumberFieldFractionalIdeal_rel
from sage.libs.all import pari, pari_gen

QQ = rational_field.RationalField()
ZZ = integer_ring.IntegerRing()
RIF = sage.rings.real_mpfi.RealIntervalField()
CIF = sage.rings.complex_interval_field.ComplexIntervalField()
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.real_lazy import RLF, CLF

# from sage.rings.number_field.number_field import is_AbsoluteNumberField
# from sage.rings.number_field.number_field import is_QuadraticField
# from sage.rings.number_field.number_field import is_CyclotomicField

def is_RelativeNumberField(x):
    """
    Return True if x is a relative number field.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field_rel import is_RelativeNumberField
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

class NumberField_relative(NumberField_generic):
    """
    EXAMPLES:
        sage: K.<a> = NumberField(x^3 - 2)
        sage: t = K['x'].gen()
        sage: L.<b> = K.extension(t^2+t+a); L
        Number Field in b with defining polynomial x^2 + x + a over its base field
    """
    def __init__(self, base, polynomial, name,
                 latex_name=None, names=None, check=True, embedding=None):
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
            <class 'sage.rings.number_field.number_field_rel.NumberField_relative'>

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
        if embedding is not None:
            raise NotImplementedError, "Embeddings not implemented for relative number fields"
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
                                     latex_name=latex_name, check=False, embedding=embedding)

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
                                          self.latex_variable_name(), self.gen_embedding())

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

    def _element_constructor_(self, x):
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

    def _coerce_map_from_(self, R):
        """
        Canonical implicit coercion of x into self.

        Elements of this field canonically coerce in, as does anything
        that coerces into the base field of this field.

        EXAMPLES:
            sage: k.<a> = NumberField([x^5 + 2, x^7 + 3])
            sage: b = k(k.base_field().gen())
            sage: b = k.coerce(k.base_field().gen())
            sage: b^7
            -3
            sage: k.coerce(2/3)
            2/3
            sage: c = a + b  # this works
        """
        if R in [int, long, ZZ, QQ, self.base_field()]:
            return self._generic_convert_map(R)
        from sage.rings.number_field.order import is_NumberFieldOrder
        if is_NumberFieldOrder(R) and R.number_field() is self:
            return self._generic_convert_map(R)
        mor = self.base_field().coerce_map_from(R)
        if mor is not None:
            return self.coerce_map_from(self.base_field()) * mor

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
        element = self.base_field().coerce(element)
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
            Isomorphism map:
              From: Vector space of dimension 9 over Rational Field
              To:   Number Field in a with defining polynomial x^3 + 3 over its base field
            sage: to_V
            Isomorphism map:
              From: Number Field in a with defining polynomial x^3 + 3 over its base field
              To:   Vector space of dimension 9 over Rational Field
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
            Isomorphism map:
              From: Number Field in c with defining polynomial x^8 + 8*x^6 + 30*x^4 - 40*x^2 + 49
              To:   Number Field in a with defining polynomial x^4 + 3 over its base field
            sage: from_L(c)
            a - b
            sage: to_L
            Isomorphism map:
              From: Number Field in a with defining polynomial x^4 + 3 over its base field
              To:   Number Field in c with defining polynomial x^8 + 8*x^6 + 30*x^4 - 40*x^2 + 49
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
        Given an element in self or an embedding of a subfield into self,
        return a relative number field $K$ isomorphic to self that is relative
        over the absolute field $\QQ(\alpha)$ or the domain of $alpha$, along
        with isomorphisms from $K$ to self and from self to K.

        INPUT:
            alpha -- an element of self, or an embedding of a subfield into self
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

            Now suppose we have K below L below M:

            sage: M = NumberField(x^8 + 2, 'a'); M
            Number Field in a with defining polynomial x^8 + 2
            sage: L, L_into_M, _ = M.subfields(4)[0]; L
            Number Field in a0 with defining polynomial x^4 + 2
            sage: K, K_into_L, _ = L.subfields(2)[0]; K
            Number Field in a00 with defining polynomial x^2 + 2
            sage: K_into_M = L_into_M * K_into_L

            sage: L_over_K = L.relativize(K_into_L, 'c'); L_over_K
            Number Field in c0 with defining polynomial x^2 + a00 over its base field
            sage: L_over_K_to_L, L_to_L_over_K = L_over_K.structure()
            sage: M_over_L_over_K = M.relativize(L_into_M * L_over_K_to_L, 'd'); M_over_L_over_K
            Number Field in d0 with defining polynomial x^2 + c0 over its base field
            sage: M_over_L_over_K.base_field() is L_over_K
            True

            Let's test relativizing a degree 6 field over its degree 2 and
            degree 3 subfields, using both an explicit element

            sage: K.<a> = NumberField(x^6 + 2); K
            Number Field in a with defining polynomial x^6 + 2
            sage: K2, K2_into_K, _ = K.subfields(2)[0]; K2
            Number Field in a0 with defining polynomial x^2 + 2
            sage: K3, K3_into_K, _ = K.subfields(3)[0]; K3
            Number Field in a0 with defining polynomial x^3 - 2

            Here we explicitly relativize over an element of K2 (not the
            generator):

            sage: L = K.relativize(K3_into_K, 'b'); L
            Number Field in b0 with defining polynomial x^2 + a0 over its base field
            sage: L_to_K, K_to_L = L.structure()
            sage: L_over_K2 = L.relativize(K_to_L(K2_into_K(K2.gen() + 1)), 'c'); L_over_K2
            Number Field in c0 with defining polynomial x^3 - c1 + 1 over its base field
            sage: L_over_K2.base_field()
            Number Field in c1 with defining polynomial x^2 - 2*x + 3

            Here we use a morphism to preserve the base field information:

            sage: K2_into_L = K_to_L * K2_into_K
            sage: L_over_K2 = L.relativize(K2_into_L, 'c'); L_over_K2
            Number Field in c0 with defining polynomial x^3 - a0 over its base field
            sage: L_over_K2.base_field() is K2
            True
        """
        K = self.absolute_field('a')
        from_K, to_K = K.structure()

        if is_Map(alpha):
            # alpha is an embedding of a subfield into self; compose to get an
            # embedding of a subfield into the absolute field
            beta = to_K * alpha
        else:
            # alpha is an element coercible into self
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

def NumberField_relative_v1(base_field, poly, name, latex_name, canonical_embedding=None):
    """
    This is used in pickling relative fields.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field_rel import NumberField_relative_v1
        sage: R.<x> = CyclotomicField(3)[]
        sage: NumberField_relative_v1(CyclotomicField(3), x^2 + 7, 'a', 'a')
        Number Field in a with defining polynomial x^2 + 7 over its base field
    """
    return NumberField_relative(base_field, poly, name, latex_name, check=False, embedding=canonical_embedding)

NumberField_extension_v1 = NumberField_relative_v1  # historical reasons only
