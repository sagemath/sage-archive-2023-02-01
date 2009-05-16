"""
Number Field Elements

AUTHORS:

- William Stein: version before it got Cython'd

- Joel B. Mohler (2007-03-09): First reimplementation in Cython

- William Stein (2007-09-04): add doctests

- Robert Bradshaw (2007-09-15): specialized classes for relative and
  absolute elements
"""

# TODO -- relative extensions need to be completely rewritten, so one
# can get easy access to representation of elements in their relative
# form.  Functions like matrix below can't be done until relative
# extensions are re-written this way.  Also there needs to be class
# hierarchy for number field elements, integers, etc.  This is a
# nontrivial project, and it needs somebody to attack it.  I'm amazed
# how long this has gone unattacked.

# Relative elements need to be a derived class or something.  This is
# terrible as it is now.

#*****************************************************************************
#       Copyright (C) 2004, 2007 William Stein <wstein@gmail.com>
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

include '../../ext/interrupt.pxi'
include '../../ext/python_int.pxi'
include "../../ext/stdsage.pxi"

import sage.rings.field_element
import sage.rings.infinity
import sage.rings.polynomial.polynomial_element
import sage.rings.rational_field
import sage.rings.rational
import sage.rings.integer_ring
import sage.rings.integer
import sage.rings.arith

import number_field

from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.rings.integer_ring cimport IntegerRing_class

from sage.modules.free_module_element import vector

from sage.libs.all import pari_gen
from sage.libs.pari.gen import PariError
from sage.structure.element cimport Element, generic_power_c

QQ = sage.rings.rational_field.QQ
ZZ = sage.rings.integer_ring.ZZ
Integer_sage = sage.rings.integer.Integer

# this is a threshold for the charpoly() methods in this file
# for degrees <= this threshold, pari is used
# for degrees > this threshold, sage matrices are used
# the value was decided by running a tuning script on a number of
# architectures; you can find this script attached to trac
# ticket 5213
TUNE_CHARPOLY_NF = 25

def is_NumberFieldElement(x):
    """
    Return True if x is of type NumberFieldElement, i.e., an element of
    a number field.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_element import is_NumberFieldElement
        sage: is_NumberFieldElement(2)
        False
        sage: k.<a> = NumberField(x^7 + 17*x + 1)
        sage: is_NumberFieldElement(a+1)
        True
    """
    return PY_TYPE_CHECK(x, NumberFieldElement)

def __create__NumberFieldElement_version0(parent, poly):
    """
    Used in unpickling elements of number fields.

    EXAMPLES:

    Since this is just used in unpickling, we unpickle.

    ::

        sage: k.<a> = NumberField(x^3 - 2)
        sage: loads(dumps(a+1)) == a + 1
        True
    """
    return NumberFieldElement(parent, poly)

def __create__NumberFieldElement_version1(parent, cls, poly):
    """
    Used in unpickling elements of number fields.

    EXAMPLES:

    Since this is just used in unpickling, we unpickle.

    ::

        sage: k.<a> = NumberField(x^3 - 2)
        sage: loads(dumps(a+1)) == a + 1
        True
    """
    return cls(parent, poly)

def _inverse_mod_generic(elt, I):
    r"""
    Return an inverse of elt modulo the given ideal. This is a separate
    function called from each of the OrderElement_xxx classes, since
    otherwise we'd have to have the same code three times over (there
    is no OrderElement_generic class - no multiple inheritance). See
    trac 4190.

    EXAMPLES::

        sage: OE = NumberField(x^3 - x + 2, 'w').ring_of_integers()
        sage: w = OE.ring_generators()[0]
        sage: from sage.rings.number_field.number_field_element import _inverse_mod_generic
        sage: _inverse_mod_generic(w, 13*OE)
        -7*w^2 - 13*w + 7
    """
    from sage.matrix.constructor import matrix
    R = elt.parent()
    n = R.absolute_degree()
    I = R.ideal_monoid()(I)
    if not I.is_integral():
        raise TypeError, "inverse modulo non-integral ideals not defined"
    m = matrix(ZZ, [R.coordinates(y) for y in I.integral_basis()] + [R.coordinates(elt*s) for s in R.gens()])
    a,b = m.echelon_form(transformation=True)
    if (a[0:n] != 1):
        raise ZeroDivisionError, "%s is not invertible modulo %s" % (elt, I)
    v = R.coordinates(1)
    y = sum([b[j,i+n] * R.gens()[i] * v[j] for i in xrange(n) for j in xrange(n)])
    return y

__pynac_pow = False

cdef class NumberFieldElement(FieldElement):
    """
    An element of a number field.

    EXAMPLES::

        sage: k.<a> = NumberField(x^3 + x + 1)
        sage: a^3
        -a - 1
    """
    cdef _new(self):
        """
        Quickly creates a new initialized NumberFieldElement with the same
        parent as self.
        """
        cdef NumberFieldElement x
        x = <NumberFieldElement>PY_NEW_SAME_TYPE(self)
        x._parent = self._parent
        x.__fld_numerator = self.__fld_numerator
        x.__fld_denominator = self.__fld_denominator
        return x

    cdef number_field(self):
        return self._parent

    def _number_field(self):
        return self.number_field()

    def __init__(self, parent, f):
        """
        INPUT:


        -  ``parent`` - a number field

        -  ``f`` - defines an element of a number field.


        EXAMPLES:

        The following examples illustrate creation of elements of
        number fields, and some basic arithmetic.

        First we define a polynomial over Q.

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^2 + 1

        Next we use f to define the number field.

        ::

            sage: K.<a> = NumberField(f); K
            Number Field in a with defining polynomial x^2 + 1
            sage: a = K.gen()
            sage: a^2
            -1
            sage: (a+1)^2
            2*a
            sage: a^2
            -1
            sage: z = K(5); 1/z
            1/5

        We create a cube root of 2.

        ::

            sage: K.<b> = NumberField(x^3 - 2)
            sage: b = K.gen()
            sage: b^3
            2
            sage: (b^2 + b + 1)^3
            12*b^2 + 15*b + 19

        This example illustrates save and load::

            sage: K.<a> = NumberField(x^17 - 2)
            sage: s = a^15 - 19*a + 3
            sage: loads(s.dumps()) == s
            True
        """
        sage.rings.field_element.FieldElement.__init__(self, parent)
        self.__fld_numerator, self.__fld_denominator = parent.absolute_polynomial_ntl()

        cdef ZZ_c coeff
        if isinstance(f, (int, long, Integer_sage)):
            # set it up and exit immediately
            # fast pathway
            (<Integer>ZZ(f))._to_ZZ(&coeff)
            ZZX_SetCoeff( self.__numerator, 0, coeff )
            ZZ_conv_from_int( self.__denominator, 1 )
            return

        elif isinstance(f, NumberFieldElement):
            if PY_TYPE(self) is PY_TYPE(f):
                self.__numerator = (<NumberFieldElement>f).__numerator
                self.__denominator = (<NumberFieldElement>f).__denominator
                return
            else:
                f = f.polynomial()

        from sage.rings.number_field import number_field_rel
        if isinstance(parent, number_field_rel.NumberField_relative):
            ppr = parent.base_field().polynomial_ring()
        else:
            ppr = parent.polynomial_ring()

        cdef long i
        if isinstance(f, pari_gen):
            if f.type() == "t_COL":
                newf = ppr(0)
                Zbasis = self.number_field().pari_nf().getattr('zk')
                # Note that this integral basis is not the same as that returned by parent.integral_basis() !
                for i from 0 <= i < parent.degree():
                    if f[i] != 0:
                        newf += QQ(f[i]) * ppr(Zbasis[i])
                f = newf
            else:
                if f.type() == "t_POLMOD":
                    f = f.lift()
                if f.type() in ["t_INT", "t_FRAC", "t_POL"]:
                    f = ppr(f)
                else:
                    raise TypeError, "Unsupported Pari type"
        f = ppr(f)
        if f.degree() >= parent.absolute_degree():
            from sage.rings.number_field import number_field_rel
            if isinstance(parent, number_field_rel.NumberField_relative):
                f %= ppr(parent.absolute_polynomial())
            else:
                f %= parent.polynomial()

        # Set Denominator
        den = f.denominator()
        (<Integer>ZZ(den))._to_ZZ(&self.__denominator)

        num = f * den
        for i from 0 <= i <= num.degree():
            (<Integer>ZZ(num[i]))._to_ZZ(&coeff)
            ZZX_SetCoeff( self.__numerator, i, coeff )

    def __new__(self, parent = None, f = None):
        ZZX_construct(&self.__numerator)
        ZZ_construct(&self.__denominator)

    def __dealloc__(self):
        ZZX_destruct(&self.__numerator)
        ZZ_destruct(&self.__denominator)

    def _lift_cyclotomic_element(self, new_parent, bint check=True, int rel=0):
        """
        Creates an element of the passed field from this field. This is
        specific to creating elements in a cyclotomic field from elements
        in another cyclotomic field, in the case that
        self.number_field()._n() divides new_parent()._n(). This
        function aims to make this common coercion extremely fast!

        More general coercion (i.e. of zeta6 into CyclotomicField(3)) is
        implemented in the _coerce_from_other_cyclotomic_field method
        of a CyclotomicField.

        EXAMPLES::

            sage: C.<zeta5>=CyclotomicField(5)
            sage: CyclotomicField(10)(zeta5+1)  # The function _lift_cyclotomic_element does the heavy lifting in the background
            zeta10^2 + 1
            sage: (zeta5+1)._lift_cyclotomic_element(CyclotomicField(10))  # There is rarely a purpose to call this function directly
            zeta10^2 + 1
            sage: cf4 = CyclotomicField(4)
            sage: cf1 = CyclotomicField(1) ; one = cf1.0
            sage: cf4(one)
            1
            sage: type(cf4(1))
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
            sage: cf33 = CyclotomicField(33) ; z33 = cf33.0
            sage: cf66 = CyclotomicField(66) ; z66 = cf66.0
            sage: z33._lift_cyclotomic_element(cf66)
            zeta66^2
            sage: z66._lift_cyclotomic_element(cf33)
            Traceback (most recent call last):
            ...
            TypeError: The zeta_order of the new field must be a multiple of the zeta_order of the original.
            sage: cf33(z66)
            -zeta33^17

        AUTHORS:

        - Joel B. Mohler

        - Craig Citro (fixed behavior for different representation of
          quadratic field elements)
        """
        if check:
            if not isinstance(self.number_field(), number_field.NumberField_cyclotomic) \
                   or not isinstance(new_parent, number_field.NumberField_cyclotomic):
                raise TypeError, "The field and the new parent field must both be cyclotomic fields."

        if rel == 0:
            small_order = self.number_field()._n()
            large_order = new_parent._n()

            try:
                rel = ZZ(large_order / small_order)
            except TypeError:
                raise TypeError, "The zeta_order of the new field must be a multiple of the zeta_order of the original."

        ## degree 2 is handled differently, because elements are
        ## represented differently
        if new_parent.degree() == 2:
            return new_parent._element_class(new_parent, self)

        cdef NumberFieldElement x = <NumberFieldElement>PY_NEW_SAME_TYPE(self)
        x._parent = <ParentWithBase>new_parent
        x.__fld_numerator, x.__fld_denominator = new_parent.polynomial_ntl()
        x.__denominator = self.__denominator
        cdef ZZX_c result
        cdef ZZ_c tmp
        cdef int i
        cdef ntl_ZZX _num
        cdef ntl_ZZ _den
        for i from 0 <= i <= ZZX_deg(self.__numerator):
            tmp = ZZX_coeff(self.__numerator, i)
            ZZX_SetCoeff(result, i*rel, tmp)
        ZZX_rem(x.__numerator, result, x.__fld_numerator.x)
        return x

    def __reduce__(self):
        """
        Used in pickling number field elements.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 17*x^2 + 1)
            sage: t = a.__reduce__(); t
            (<built-in function __create__NumberFieldElement_version1>, (Number Field in a with defining polynomial x^3 - 17*x^2 + 1, <type 'sage.rings.number_field.number_field_element.NumberFieldElement_absolute'>, x))
            sage: t[0](*t[1]) == a
            True
        """
        return __create__NumberFieldElement_version1, \
               (self.number_field(), type(self), self.polynomial())

    def __repr__(self):
        """
        String representation of this number field element, which is just a
        polynomial in the generator.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 2)
            sage: b = (2/3)*a + 3/5
            sage: b.__repr__()
            '2/3*a + 3/5'
        """
        x = self.polynomial()
        K = self.number_field()
        return str(x).replace(x.parent().variable_name(), K.variable_name())

    def _im_gens_(self, codomain, im_gens):
        """
        This is used in computing homomorphisms between number fields.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 - 2)
            sage: m.<b> = NumberField(x^4 - 2)
            sage: phi = k.hom([b^2])
            sage: phi(a+1)
            b^2 + 1
            sage: (a+1)._im_gens_(m, [b^2])
            b^2 + 1
        """
        # NOTE -- if you ever want to change this so relative number
        # fields are in terms of a root of a poly.  The issue is that
        # elements of a relative number field are represented in terms
        # of a generator for the absolute field.  However the morphism
        # gives the image of gen, which need not be a generator for
        # the absolute field.  The morphism has to be *over* the
        # relative element.
        return codomain(self.polynomial()(im_gens[0]))

    def _latex_(self):
        """
        Returns the latex representation for this element.

        EXAMPLES::

            sage: C,zeta12=CyclotomicField(12).objgen()
            sage: latex(zeta12^4-zeta12)
            \zeta_{12}^{2} - \zeta_{12} - 1
        """
        return self.polynomial()._latex_(name=self.number_field().latex_variable_name())

    def _gap_init_(self):
        """
        Return gap string representation of self.

        EXAMPLES::

            sage: F=CyclotomicField(8)
            sage: p=F.gen()^2+2*F.gen()-3
            sage: p
            zeta8^2 + 2*zeta8 - 3
            sage: p._gap_init_() # The variable name $sage2 belongs to the gap(F) and is somehow random
            'GeneratorsOfField($sage2)[1]^2 + 2*GeneratorsOfField($sage2)[1] - 3'
            sage: gap(p._gap_init_())
            (-3+2*zeta8+zeta8^2)
        """
        s = self.__repr__()
        return s.replace(str(self.parent().gen()), 'GeneratorsOfField(%s)[1]'%sage.interfaces.gap.gap(self.parent()).name())


    def _pari_(self, var='x'):
        raise NotImplementedError, "NumberFieldElement sub-classes must override _pari_"

    def _pari_init_(self, var='x'):
        """
        Return GP/PARI string representation of self. This is used for
        converting this number field element to GP/PARI. The returned
        string defines a pari Mod in the variable is var, which is by
        default 'x' - not the name of the generator of the number field.

        INPUT:


        -  ``var`` - (default: 'x') the variable of the pari
           Mod.


        EXAMPLES::

            sage: K.<a> = NumberField(x^5 - x - 1)
            sage: ((1 + 1/3*a)^4)._pari_init_()
            'Mod(1/81*x^4 + 4/27*x^3 + 2/3*x^2 + 4/3*x + 1, x^5 - x - 1)'
            sage: ((1 + 1/3*a)^4)._pari_init_('a')
            'Mod(1/81*a^4 + 4/27*a^3 + 2/3*a^2 + 4/3*a + 1, a^5 - a - 1)'

        Note that _pari_init_ can fail because of reserved words in
        PARI, and since it actually works by obtaining the PARI
        representation of something.

        ::

            sage: K.<theta> = NumberField(x^5 - x - 1)
            sage: b = (1/2 - 2/3*theta)^3; b
            -8/27*theta^3 + 2/3*theta^2 - 1/2*theta + 1/8
            sage: b._pari_init_('theta')
            Traceback (most recent call last):
            ...
            PariError: unexpected character (2)

        Fortunately pari_init returns everything in terms of x by
        default.

        ::

            sage: pari(b)
            Mod(-8/27*x^3 + 2/3*x^2 - 1/2*x + 1/8, x^5 - x - 1)
        """
        return repr(self._pari_(var=var))
##         if var == None:
##             var = self.parent().variable_name()
##         if isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_relative):
##             f = self.polynomial()._pari_()
##             g = str(self.parent().pari_relative_polynomial())
##             base = self.parent().base_ring()
##             gsub = base.gen()._pari_()
##             gsub = str(gsub).replace(base.variable_name(), "y")
##             g = g.replace("y", gsub)
##         else:
##             f = str(self.polynomial()).replace("x",var)
##             g = str(self.parent().polynomial()).replace("x",var)
##         return 'Mod(%s, %s)'%(f,g)

    def __getitem__(self, n):
        """
        Return the n-th coefficient of this number field element, written
        as a polynomial in the generator.

        Note that `n` must be between 0 and `d-1`, where
        `d` is the degree of the number field.

        EXAMPLES::

            sage: m.<b> = NumberField(x^4 - 1789)
            sage: c = (2/3-4/5*b)^3; c
            -64/125*b^3 + 32/25*b^2 - 16/15*b + 8/27
            sage: c[0]
            8/27
            sage: c[2]
            32/25
            sage: c[3]
            -64/125

        We illustrate bounds checking::

            sage: c[-1]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1.
            sage: c[4]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1.

        The list method implicitly calls ``__getitem__``::

            sage: list(c)
            [8/27, -16/15, 32/25, -64/125]
            sage: m(list(c)) == c
            True
        """
        if n < 0 or n >= self.number_field().degree():     # make this faster.
            raise IndexError, "index must be between 0 and degree minus 1."
        return self.polynomial()[n]

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef NumberFieldElement _right = right
        return not (ZZX_equal(left.__numerator, _right.__numerator) and ZZ_equal(left.__denominator, _right.__denominator))

    def __abs__(self):
        r"""
        Return the numerical absolute value of this number field element
        with respect to the first archimedean embedding, to double
        precision.

        This is the ``abs( )`` Python function. If you want a
        different embedding or precision, use
        ``self.abs(...)``.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 2)
            sage: abs(a)
            1.25992104989
            sage: abs(a)^3
            2.0
            sage: a.abs(prec=128)
            1.2599210498948731647672106072782283506
        """
        return self.abs(prec=53, i=0)

    def abs(self, prec=53, i=0):
        """
        Return the absolute value of this element with respect to the ith
        complex embedding of parent, to the given precision.

        If prec is 53 (the default), then the complex double field is used;
        otherwise the arbitrary precision (but slow) complex field is
        used.

        INPUT:


        -  ``prec`` - (default: 53) integer bits of precision

        -  ``i`` - (default: ) integer, which embedding to
           use


        EXAMPLES::

            sage: z = CyclotomicField(7).gen()
            sage: abs(z)
            1.0
            sage: abs(z^2 + 17*z - 3)
            16.06044268
            sage: K.<a> = NumberField(x^3+17)
            sage: abs(a)
            2.57128159066
            sage: a.abs(prec=100)
            2.5712815906582353554531872087
            sage: a.abs(prec=100,i=1)
            2.5712815906582353554531872087
            sage: a.abs(100, 2)
            2.5712815906582353554531872087

        Here's one where the absolute value depends on the embedding.

        ::

            sage: K.<b> = NumberField(x^2-2)
            sage: a = 1 + b
            sage: a.abs(i=0)
            0.414213562373
            sage: a.abs(i=1)
            2.41421356237
        """
        P = self.number_field().complex_embeddings(prec)[i]
        return abs(P(self))

    def coordinates_in_terms_of_powers(self):
        r"""
        Let `\alpha` be self. Return a Python function that takes
        any element of the parent of self in `\QQ(\alpha)`
        and writes it in terms of the powers of `\alpha`:
        `1, \alpha, \alpha^2, ...`.

        (NOT CACHED).

        EXAMPLES:

        This function allows us to write elements of a number
        field in terms of a different generator without having to construct
        a whole separate number field.

        ::

            sage: y = polygen(QQ,'y'); K.<beta> = NumberField(y^3 - 2); K
            Number Field in beta with defining polynomial y^3 - 2
            sage: alpha = beta^2 + beta + 1
            sage: c = alpha.coordinates_in_terms_of_powers(); c
            Coordinate function that writes elements in terms of the powers of beta^2 + beta + 1
            sage: c(beta)
            [-2, -3, 1]
            sage: c(alpha)
            [0, 1, 0]
            sage: c((1+beta)^5)
            [3, 3, 3]
            sage: c((1+beta)^10)
            [54, 162, 189]

        This function works even if self only generates a subfield of this
        number field.

        ::

            sage: k.<a> = NumberField(x^6 - 5)
            sage: alpha = a^3
            sage: c = alpha.coordinates_in_terms_of_powers()
            sage: c((2/3)*a^3 - 5/3)
            [-5/3, 2/3]
            sage: c
            Coordinate function that writes elements in terms of the powers of a^3
            sage: c(a)
            Traceback (most recent call last):
            ...
            ArithmeticError: vector is not in free module
        """
        K = self.number_field()
        V, from_V, to_V = K.absolute_vector_space()
        h = K(1)
        B = [to_V(h)]
        f = self.minpoly()
        for i in range(f.degree()-1):
            h *= self
            B.append(to_V(h))
        W = V.span_of_basis(B)
        return CoordinateFunction(self, W, to_V)

    def complex_embeddings(self, prec=53):
        """
        Return the images of this element in the floating point complex
        numbers, to the given bits of precision.

        INPUT:


        -  ``prec`` - integer (default: 53) bits of precision


        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 2)
            sage: a.complex_embeddings()
            [-0.629960524947 - 1.09112363597*I, -0.629960524947 + 1.09112363597*I, 1.25992104989]
            sage: a.complex_embeddings(10)
            [-0.63 - 1.1*I, -0.63 + 1.1*I, 1.3]
            sage: a.complex_embeddings(100)
            [-0.62996052494743658238360530364 - 1.0911236359717214035600726142*I, -0.62996052494743658238360530364 + 1.0911236359717214035600726142*I, 1.2599210498948731647672106073]
        """
        phi = self.number_field().complex_embeddings(prec)
        return [f(self) for f in phi]

    def complex_embedding(self, prec=53, i=0):
        """
        Return the i-th embedding of self in the complex numbers, to the
        given precision.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 2)
            sage: a.complex_embedding()
            -0.629960524947 - 1.09112363597*I
            sage: a.complex_embedding(10)
            -0.63 - 1.1*I
            sage: a.complex_embedding(100)
            -0.62996052494743658238360530364 - 1.0911236359717214035600726142*I
            sage: a.complex_embedding(20, 1)
            -0.62996 + 1.0911*I
            sage: a.complex_embedding(20, 2)
            1.2599
        """
        return self.number_field().complex_embeddings(prec)[i](self)

    def _complex_double_(self, CDF):
        """
        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: abs(CDF(a))
	    1.0
        """
        return CDF(self.complex_embedding())

    def is_totally_positive(self):
        """
        Returns True if self is positive for all real embeddings of its
        parent number field. We do nothing at complex places, so e.g. any
        element of a totally complex number field will return True.

        EXAMPLES::

            sage: F.<b> = NumberField(x^3-3*x-1)
            sage: b.is_totally_positive()
            False
            sage: (b^2).is_totally_positive()
            True
        """
        for v in self.number_field().real_embeddings():
            if v(self) <= 0:
                return False
        return True

    def is_square(self, root=False):
        """
        Return True if self is a square in its parent number field and
        otherwise return False.

        INPUT:


        -  ``root`` - if True, also return a square root (or
           None if self is not a perfect square)


        EXAMPLES::

            sage: m.<b> = NumberField(x^4 - 1789)
            sage: b.is_square()
            False
            sage: c = (2/3*b + 5)^2; c
            4/9*b^2 + 20/3*b + 25
            sage: c.is_square()
            True
            sage: c.is_square(True)
            (True, 2/3*b + 5)

        We also test the functional notation.

        ::

            sage: is_square(c, True)
            (True, 2/3*b + 5)
            sage: is_square(c)
            True
            sage: is_square(c+1)
            False
        """
        v = self.sqrt(all=True)
        t = len(v) > 0
        if root:
            if t:
                return t, v[0]
            else:
                return False, None
        else:
            return t

    def sqrt(self, all=False):
        """
        Returns the square root of this number in the given number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 3)
            sage: K(3).sqrt()
            a
            sage: K(3).sqrt(all=True)
            [a, -a]
            sage: K(a^10).sqrt()
            9*a
            sage: K(49).sqrt()
            7
            sage: K(1+a).sqrt()
            Traceback (most recent call last):
            ...
            ValueError: a + 1 not a square in Number Field in a with defining polynomial x^2 - 3
            sage: K(0).sqrt()
            0
            sage: K((7+a)^2).sqrt(all=True)
            [a + 7, -a - 7]

        ::

            sage: K.<a> = CyclotomicField(7)
            sage: a.sqrt()
            a^4

        ::

            sage: K.<a> = NumberField(x^5 - x + 1)
            sage: (a^4 + a^2 - 3*a + 2).sqrt()
            a^3 - a^2

        ALGORITHM: Use Pari to factor `x^2` - ``self``
        in K.
        """
        # For now, use pari's factoring abilities
        R = self.number_field()['t']
        f = R([-self, 0, 1])
        roots = f.roots()
        if all:
            return [r[0] for r in roots]
        elif len(roots) > 0:
            return roots[0][0]
        else:
            try:
                # This is what integers, rationals do...
                from sage.all import SR, sqrt
                return sqrt(SR(self))
            except TypeError:
                raise ValueError, "%s not a square in %s"%(self, self._parent)

    def nth_root(self, n, all=False):
        r"""
        Return an nth root of self in the given number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4-7)
            sage: K(7).nth_root(2)
            a^2
            sage: K((a-3)^5).nth_root(5)
            a - 3

        ALGORITHM: Use Pari to factor `x^n` - ``self``
        in K.
        """
        R = self.number_field()['t']
        if not self:
            return [self] if all else self
        f = (R.gen(0) << (n-1)) - self
        roots = f.roots()
        if all:
            return [r[0] for r in roots]
        elif len(roots) > 0:
            return roots[0][0]
        else:
            raise ValueError, "%s not a %s-th root in %s"%(self, n, self._parent)

    def __pow__(base, exp, dummy):
        """
        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: sqrt2^2
            2
            sage: sqrt2^5
            4*sqrt2
            sage: (1+sqrt2)^100
            66992092050551637663438906713182313772*sqrt2 + 94741125149636933417873079920900017937
            sage: (1+sqrt2)^-1
            sqrt2 - 1

        If the exponent is not integral, perform this operation in
        the symbolic ring::

            sage: sqrt2^(1/5)
            2^(1/10)
            sage: sqrt2^sqrt2
            2^(1/2*sqrt(2))

        TESTS::

            sage: 2^I
            2^I
        """
        if (PY_TYPE_CHECK(base, NumberFieldElement) and
            (PY_TYPE_CHECK(exp, Integer) or PY_TYPE_CHECK_EXACT(exp, int) or exp in ZZ)):
            return generic_power_c(base, exp, None)
        else:
            from sage.symbolic.power_helper import try_symbolic_power
            return try_symbolic_power(base, exp)

    cdef void _reduce_c_(self):
        """
        Pull out common factors from the numerator and denominator!
        """
        cdef ZZ_c gcd
        cdef ZZ_c t1
        cdef ZZX_c t2
        ZZX_content(t1, self.__numerator)
        ZZ_GCD(gcd, t1, self.__denominator)
        if ZZ_sign(gcd) != ZZ_sign(self.__denominator):
            ZZ_negate(t1, gcd)
            gcd = t1
        ZZX_div_ZZ(t2, self.__numerator, gcd)
        ZZ_div(t1, self.__denominator, gcd)
        self.__numerator = t2
        self.__denominator = t1

    cpdef ModuleElement _add_(self, ModuleElement right):
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        x = self._new()
        ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
        cdef ZZX_c t1, t2
        ZZX_mul_ZZ(t1, self.__numerator, _right.__denominator)
        ZZX_mul_ZZ(t2, _right.__numerator, self.__denominator)
        ZZX_add(x.__numerator, t1, t2)
        x._reduce_c_()
        return x

    cpdef ModuleElement _sub_(self, ModuleElement right):
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        x = self._new()
        ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
        cdef ZZX_c t1, t2
        ZZX_mul_ZZ(t1, self.__numerator, _right.__denominator)
        ZZX_mul_ZZ(t2, _right.__numerator, self.__denominator)
        ZZX_sub(x.__numerator, t1, t2)
        x._reduce_c_()
        return x

    cpdef RingElement _mul_(self, RingElement right):
        """
        Returns the product of self and other as elements of a number
        field.

        EXAMPLES::

            sage: C.<zeta12>=CyclotomicField(12)
            sage: zeta12*zeta12^11
            1
            sage: G.<a> = NumberField(x^3 + 2/3*x + 1)
            sage: a^3
            -2/3*a - 1
            sage: a^3+a
            1/3*a - 1
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        cdef ZZX_c temp
        cdef ZZ_c temp1
        x = self._new()
        _sig_on
        # MulMod doesn't handle non-monic polynomials.
        # Therefore, we handle the non-monic case entirely separately.

        if ZZ_IsOne(ZZX_LeadCoeff(self.__fld_numerator.x)):
            ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
            ZZX_MulMod(x.__numerator, self.__numerator, _right.__numerator, self.__fld_numerator.x)
        else:
            ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
            ZZX_mul(x.__numerator, self.__numerator, _right.__numerator)
            if ZZX_deg(x.__numerator) >= ZZX_deg(self.__fld_numerator.x):
                ZZX_mul_ZZ( x.__numerator, x.__numerator, self.__fld_denominator.x )
                ZZX_mul_ZZ( temp, self.__fld_numerator.x, x.__denominator )
                ZZ_power(temp1,ZZX_LeadCoeff(temp),ZZX_deg(x.__numerator)-ZZX_deg(self.__fld_numerator.x)+1)
                ZZX_PseudoRem(x.__numerator, x.__numerator, temp)
                ZZ_mul(x.__denominator, x.__denominator, self.__fld_denominator.x)
                ZZ_mul(x.__denominator, x.__denominator, temp1)
        _sig_off
        x._reduce_c_()
        return x

        #NOTES: In LiDIA, they build a multiplication table for the
        #number field, so it's not necessary to reduce modulo the
        #defining polynomial every time:
        #     src/number_fields/algebraic_num/order.cc: compute_table
        # but asymptotically fast poly multiplication means it's
        # actually faster to *not* build a table!?!

    cpdef RingElement _div_(self, RingElement right):
        """
        Returns the quotient of self and other as elements of a number
        field.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: 1/I
            -I
            sage: I/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        ::

            sage: G.<a> = NumberField(x^3 + 2/3*x + 1)
            sage: a/a
            1
            sage: 1/a
            -a^2 - 2/3
            sage: a/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Number field element division by zero
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        cdef ZZX_c inv_num
        cdef ZZ_c inv_den
        cdef ZZX_c temp
        cdef ZZ_c temp1
        if not _right:
            raise ZeroDivisionError, "Number field element division by zero"
        x = self._new()
        _sig_on
        _right._invert_c_(&inv_num, &inv_den)
        if ZZ_IsOne(ZZX_LeadCoeff(self.__fld_numerator.x)):
            ZZ_mul(x.__denominator, self.__denominator, inv_den)
            ZZX_MulMod(x.__numerator, self.__numerator, inv_num, self.__fld_numerator.x)
        else:
            ZZ_mul(x.__denominator, self.__denominator, inv_den)
            ZZX_mul(x.__numerator, self.__numerator, inv_num)
            if ZZX_deg(x.__numerator) >= ZZX_deg(self.__fld_numerator.x):
                ZZX_mul_ZZ( x.__numerator, x.__numerator, self.__fld_denominator.x )
                ZZX_mul_ZZ( temp, self.__fld_numerator.x, x.__denominator )
                ZZ_power(temp1,ZZX_LeadCoeff(temp),ZZX_deg(x.__numerator)-ZZX_deg(self.__fld_numerator.x)+1)
                ZZX_PseudoRem(x.__numerator, x.__numerator, temp)
                ZZ_mul(x.__denominator, x.__denominator, self.__fld_denominator.x)
                ZZ_mul(x.__denominator, x.__denominator, temp1)
        x._reduce_c_()
        _sig_off
        return x

    def __floordiv__(self, other):
        """
        Return the quotient of self and other. Since these are field
        elements the floor division is exactly the same as usual division.

        EXAMPLES::

            sage: m.<b> = NumberField(x^4 + x^2 + 2/3)
            sage: c = (1+b) // (1-b); c
            3/4*b^3 + 3/4*b^2 + 3/2*b + 1/2
            sage: (1+b) / (1-b) == c
            True
            sage: c * (1-b)
            b + 1
        """
        return self / other

    def __nonzero__(self):
        """
        Return True if this number field element is nonzero.

        EXAMPLES::

            sage: m.<b> = CyclotomicField(17)
            sage: m(0).__nonzero__()
            False
            sage: b.__nonzero__()
            True

        Nonzero is used by the bool command::

            sage: bool(b + 1)
            True
            sage: bool(m(0))
            False
        """
        return not IsZero_ZZX(self.__numerator)

    cpdef ModuleElement _neg_(self):
        cdef NumberFieldElement x
        x = self._new()
        ZZX_mul_long(x.__numerator, self.__numerator, -1)
        x.__denominator = self.__denominator
        return x

    def __copy__(self):
        cdef NumberFieldElement x
        x = self._new()
        x.__numerator = self.__numerator
        x.__denominator = self.__denominator
        return x

    def __int__(self):
        """
        Attempt to convert this number field element to a Python integer,
        if possible.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: int(1/I)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
            sage: int(I*I)
            -1

        ::

            sage: K.<a> = NumberField(x^10 - x - 1)
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
            sage: int(K(9390283))
            9390283

        The semantics are like in Python, so the value does not have to
        preserved.

        ::

            sage: int(K(393/29))
            13
        """
        return int(self.polynomial())

    def __long__(self):
        """
        Attempt to convert this number field element to a Python long, if
        possible.

        EXAMPLES::

            sage: K.<a> = NumberField(x^10 - x - 1)
            sage: long(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to long
            sage: long(K(1234))
            1234L

        The value does not have to be preserved, in the case of fractions.

        ::

            sage: long(K(393/29))
            13L
        """
        return long(self.polynomial())

    cdef void _parent_poly_c_(self, ZZX_c *num, ZZ_c *den):
        """
        I believe this function should be removed since I've put the
        pointer __fld_numerator and __fld_denominator in the element
        class. I'm not going to remove it quite yet, but feel free to
        remove it if you agree with me that it should go.
        """
        raise NotImplementedError, "NumberFieldElement subclasses must override _parent_poly_c_()"

    cdef void _invert_c_(self, ZZX_c *num, ZZ_c *den):
        """
        Computes the numerator and denominator of the multiplicative
        inverse of this element.

        Suppose that this element is x/d and the parent mod'ding polynomial
        is M/D. The NTL function XGCD( r, s, t, a, b ) computes r,s,t such
        that `r=s*a+t*b`. We compute XGCD( r, s, t, x\*D, M\*d )
        and set num=s\*D\*d den=r

        EXAMPLES:

        I'd love to, but since we are dealing with c-types, I
        can't at this level. Check __invert__ for doc-tests that rely
        on this functionality.
        """
        cdef ZZX_c t # unneeded except to be there
        cdef ZZX_c a, b
        ZZX_mul_ZZ( a, self.__numerator, self.__fld_denominator.x )
        ZZX_mul_ZZ( b, self.__fld_numerator.x, self.__denominator )
        ZZX_XGCD( den[0], num[0],  t, a, b, 1 )
        ZZX_mul_ZZ( num[0], num[0], self.__fld_denominator.x )
        ZZX_mul_ZZ( num[0], num[0], self.__denominator )

    def __invert__(self):
        """
        Returns the multiplicative inverse of self in the number field.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: ~I
            -I
            sage: (2*I).__invert__()
            -1/2*I
        """
        if IsZero_ZZX(self.__numerator):
            raise ZeroDivisionError
        cdef NumberFieldElement x
        x = self._new()
        self._invert_c_(&x.__numerator, &x.__denominator)
        x._reduce_c_()
        return x
#        K = self.parent()
#        quotient = K(1)._pari_('x') / self._pari_('x')
#        if isinstance(K, sage.rings.number_field.number_field.NumberField_relative):
#            return K(K.pari_rnf().rnfeltreltoabs(quotient))
#        else:
#            return K(quotient)

    def _integer_(self, Z=None):
        """
        Returns an integer if this element is actually an integer.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: (~I)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce -I to an integer
            sage: (2*I*I)._integer_()
            -2
        """
        if ZZX_deg(self.__numerator) >= 1:
            raise TypeError, "Unable to coerce %s to an integer"%self
        return ZZ(self._rational_())

    def _rational_(self):
        """
        Returns a rational number if this element is actually a rational
        number.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: (~I)._rational_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce -I to a rational
            sage: (I*I/2)._rational_()
            -1/2
        """
        if ZZX_deg(self.__numerator) >= 1:
            raise TypeError, "Unable to coerce %s to a rational"%self
        cdef Integer num
        num = PY_NEW(Integer)
        ZZX_getitem_as_mpz(&num.value, &self.__numerator, 0)
        return num / (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)

    def _symbolic_(self, SR):
        """
        If an embedding into CC is specified, then a representation of this
        element can be made in the symbolic ring (assuming roots of the
        minimal polynomial can be found symbolically).

        EXAMPLES::

            sage: K.<a> = QuadraticField(2)
            sage: SR(a)
            sqrt(2)
            sage: SR(3*a-5)
            3*sqrt(2) - 5
            sage: K.<a> = QuadraticField(2, embedding=-1.4)
            sage: SR(a)
            -sqrt(2)
            sage: K.<a> = NumberField(x^2 - 2)
            sage: SR(a)
            Traceback (most recent call last):
            ...
            TypeError: An embedding into RR or CC must be specified.

        Now a more complicated example::

            sage: K.<a> = NumberField(x^3 + x - 1, embedding=0.68)
            sage: b = SR(a); b
            (1/18*sqrt(3)*sqrt(31) + 1/2)^(1/3) - 1/3/(1/18*sqrt(3)*sqrt(31) + 1/2)^(1/3)

            sage: (b^3 + b - 1).simplify_radical()
            0

        Make sure we got the right one::

            sage: CC(a)
            0.682327803828019
            sage: CC(b)
            0.682327803828019

        Special case for cyclotomic fields::

            sage: K.<zeta> = CyclotomicField(19)
            sage: SR(zeta)
            e^(2/19*I*pi)
            sage: CC(zeta)
            0.945817241700635 + 0.324699469204683*I
            sage: CC(SR(zeta))
            0.945817241700635 + 0.324699469204683*I

            sage: SR(zeta^5 + 2)
            e^(10/19*I*pi) + 2

        For degree greater than 5, sometimes Galois theory prevents a
        closed-form solution.  In this case, a numerical approximation
        is used::

            sage: K.<a> = NumberField(x^5-x+1, embedding=-1)
            sage: SR(a)
            -1.1673040153

        ::

            sage: K.<a> = NumberField(x^6-x^3-1, embedding=1)
            sage: SR(a)
            1/2*(sqrt(5) + 1)^(1/3)*2^(2/3)
        """
        if self.__symbolic is None:

            K = self._parent.fraction_field()

            gen = K.gen()
            if not self is gen:
                try:
                    # share the hard work...
                    gen_image = gen._symbolic_(SR)
                    self.__symbolic = self.polynomial()(gen_image)
                    return self.__symbolic
                except TypeError:
                    pass # we may still be able to do this particular element...

            embedding = K.specified_complex_embedding()
            if embedding is None:
                raise TypeError, "An embedding into RR or CC must be specified."

            if isinstance(K, number_field.NumberField_cyclotomic):
                # solution by radicals may be difficult, but we have a closed form
                from sage.all import exp, I, pi, ComplexField, RR
                CC = ComplexField(53)
                two_pi_i = 2 * pi * I
                k = ( K._n()*CC(K.gen()).log() / CC(two_pi_i) ).real().round() # n ln z / (2 pi i)
                gen_image = exp(k*two_pi_i/K._n())
                if self is gen:
                    self.__symbolic = gen_image
                else:
                    self.__symbolic = self.polynomial()(gen_image)
            else:
                # try to solve the minpoly and choose the closest root
                poly = self.minpoly()
                roots = []
                var = SR(poly.variable_name())
                for soln in SR(poly).solve(var):
                    if soln.lhs() == var:
                        roots.append(soln.rhs())
                if len(roots) != poly.degree():
                    raise TypeError, "Unable to solve by radicals."
                from number_field_morphisms import matching_root
                from sage.rings.complex_field import ComplexField
                gen_image = matching_root(roots, self, ambient_field=ComplexField(53), margin=2)
                if gen_image is not None:
                    self.__symbolic = gen_image
                else:
                    # should be rare, e.g. if there is insufficient precision
                    raise TypeError, "Unable to determine which root in SR is this element."

        return self.__symbolic

    def galois_conjugates(self, K):
        r"""
        Return all Gal(Qbar/Q)-conjugates of this number field element in
        the field K.

        EXAMPLES:

        In the first example the conjugates are obvious::

            sage: K.<a> = NumberField(x^2 - 2)
            sage: a.galois_conjugates(K)
            [a, -a]
            sage: K(3).galois_conjugates(K)
            [3]

        In this example the field is not Galois, so we have to pass to an
        extension to obtain the Galois conjugates.

        ::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: c = a.galois_conjugates(K); c
            [a]
            sage: K.<a> = NumberField(x^3 - 2)
            sage: c = a.galois_conjugates(K.galois_closure('a1')); c
            [1/84*a1^4 + 13/42*a1, -1/252*a1^4 - 55/126*a1, -1/126*a1^4 + 8/63*a1]
            sage: c[0]^3
            2
            sage: parent(c[0])
            Number Field in a1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: parent(c[0]).is_galois()
            True

        There is only one Galois conjugate of `\sqrt[3]{2}` in
        `\QQ(\sqrt[3]{2})`.

        ::

            sage: a.galois_conjugates(K)
            [a]

        Galois conjugates of `\sqrt[3]{2}` in the field
        `\QQ(\zeta_3,\sqrt[3]{2})`::

            sage: L.<a> = CyclotomicField(3).extension(x^3 - 2)
            sage: a.galois_conjugates(L)
            [a, (-zeta3 - 1)*a, zeta3*a]
        """
        f = self.absolute_minpoly()
        g = K['x'](f)
        return [a for a,_ in g.roots()]

    def conjugate(self):
        """
        Return the complex conjugate of the number field element.
        Currently, this is implemented for cyclotomic fields and quadratic
        extensions of Q. It seems likely that there are other number fields
        for which the idea of a conjugate would be easy to compute.

        EXAMPLES::

            sage: k.<I> = QuadraticField(-1)
            sage: I.conjugate()
            -I
            sage: (I/(1+I)).conjugate()
            -1/2*I + 1/2
            sage: z6=CyclotomicField(6).gen(0)
            sage: (2*z6).conjugate()
            -2*zeta6 + 2
            sage: K.<j,b> = QQ[sqrt(-1), sqrt(2)]
            sage: j.conjugate()
            Traceback (most recent call last):
            ...
            NotImplementedError: complex conjugation is not implemented (or doesn't make sense).

        ::

            sage: K.<b> = NumberField(x^3 - 2)
            sage: b.conjugate()
            Traceback (most recent call last):
            ...
            NotImplementedError: complex conjugation is not implemented (or doesn't make sense).
        """
        coeffs = self.number_field().absolute_polynomial().list()
        if len(coeffs) == 3 and coeffs[2] == 1 and coeffs[1] == 0:
            # polynomial looks like x^2+d
            # i.e. we live in a quadratic extension of QQ
            if coeffs[0] > 0:
                gen = self.number_field().gen()
                return self.polynomial()(-gen)
            else:
                return self
        elif isinstance(self.number_field(), number_field.NumberField_cyclotomic):
            # We are in a cyclotomic field
            # Replace the generator zeta_n with (zeta_n)^(n-1)
            gen = self.number_field().gen()
            return self.polynomial()(gen ** (gen.multiplicative_order()-1))
        else:
            raise NotImplementedError, "complex conjugation is not implemented (or doesn't make sense)."

    def polynomial(self, var='x'):
        """
        Return the underlying polynomial corresponding to this number field
        element.

        The resulting polynomial is currently *not* cached.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 - x - 1)
            sage: f = (-2/3 + 1/3*a)^4; f
            1/81*a^4 - 8/81*a^3 + 8/27*a^2 - 32/81*a + 16/81
            sage: g = f.polynomial(); g
            1/81*x^4 - 8/81*x^3 + 8/27*x^2 - 32/81*x + 16/81
            sage: parent(g)
            Univariate Polynomial Ring in x over Rational Field

        Note that the result of this function is not cached (should this be
        changed?)::

            sage: g is f.polynomial()
            False
        """
        return QQ[var](self._coefficients())

    def __hash__(self):
        """
        Return hash of this number field element, which is just the
        hash of the underlying polynomial.
        """
        return hash(self.polynomial())

    def _coefficients(self):
        """
        Return the coefficients of the underlying polynomial corresponding
        to this number field element.

        OUTPUT:

        - a list whose length corresponding to the degree of this
          element written in terms of a generator.

        EXAMPLES:
        """
        coeffs = []
        cdef Integer den = (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)
        cdef Integer numCoeff
        cdef int i
        for i from 0 <= i <= ZZX_deg(self.__numerator):
            numCoeff = PY_NEW(Integer)
            ZZX_getitem_as_mpz(&numCoeff.value, &self.__numerator, i)
            coeffs.append( numCoeff / den )
        return coeffs

    cdef void _ntl_coeff_as_mpz(self, mpz_t* z, long i):
        if i > ZZX_deg(self.__numerator):
            mpz_set_ui(z[0], 0)
        else:
            ZZX_getitem_as_mpz(z, &self.__numerator, i)

    cdef void _ntl_denom_as_mpz(self, mpz_t* z):
        cdef Integer denom = (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)
        mpz_set(z[0], denom.value)

    def denominator(self):
        """
        Return the denominator of this element, which is by definition the
        denominator of the corresponding polynomial representation. I.e.,
        elements of number fields are represented as a polynomial (in
        reduced form) modulo the modulus of the number field, and the
        denominator is the denominator of this polynomial.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: a = 1/3 + (1/5)*z
            sage: print a.denominator()
            15
        """
        return (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)

    def _set_multiplicative_order(self, n):
        """
        Set the multiplicative order of this number field element.

        .. warning::

           Use with caution - only for internal use! End users should
           never call this unless they have a very good reason to do
           so.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: a._set_multiplicative_order(3)
            sage: a.multiplicative_order()
            3

        You can be evil with this so be careful. That's why the function
        name begins with an underscore.

        ::

            sage: a._set_multiplicative_order(389)
            sage: a.multiplicative_order()
            389
        """
        self.__multiplicative_order = n

    def multiplicative_order(self):
        """
        Return the multiplicative order of this number field element.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(5)
            sage: z.multiplicative_order()
            5
            sage: (-z).multiplicative_order()
            10
            sage: (1+z).multiplicative_order()
            +Infinity

            sage: x = polygen(QQ)
            sage: K.<a>=NumberField(x^40 - x^20 + 4)
            sage: u = 1/4*a^30 + 1/4*a^10 + 1/2
            sage: u.multiplicative_order()
            6
            sage: a.multiplicative_order()
            +Infinity

        An example in a relative extension:

            sage: K.<a, b> = NumberField([x^2 + x + 1, x^2 - 3])
            sage: z = (a - 1)*b/3
            sage: z.multiplicative_order()
            12
            sage: z^12==1 and z^6!=1 and z^4!=1
            True

        """
        if self.__multiplicative_order is not None:
            return self.__multiplicative_order

        one = self.number_field().one_element()
        infinity = sage.rings.infinity.infinity

        if self == one:
            self.__multiplicative_order = ZZ(1)
            return self.__multiplicative_order
        if self == -one:
            self.__multiplicative_order = ZZ(2)
            return self.__multiplicative_order

        if isinstance(self.number_field(), number_field.NumberField_cyclotomic):
            t = self.number_field()._multiplicative_order_table()
            f = self.polynomial()
            if t.has_key(f):
                self.__multiplicative_order = t[f]
                return self.__multiplicative_order
            else:
                self.__multiplicative_order = sage.rings.infinity.infinity
                return self.__multiplicative_order

        if self.is_rational_c() or not self.is_integral() or not self.norm() ==1:
            self.__multiplicative_order = infinity
            return self.__multiplicative_order

        # Now we have a unit of norm 1, and check if it is a root of unity

        n = self.number_field().zeta_order()
        if not self**n ==1:
            self.__multiplicative_order = infinity
            return self.__multiplicative_order
        from sage.groups.generic import order_from_multiple
        self.__multiplicative_order = order_from_multiple(self,n,operation='*')
        return self.__multiplicative_order

    def additive_order(self):
        r"""
        Return the additive order of this element (i.e. infinity if
        self != 0, 1 if self == 0)

        EXAMPLES::

            sage: K.<u> = NumberField(x^4 - 3*x^2 + 3)
            sage: u.additive_order()
            +Infinity
            sage: K(0).additive_order()
            1
            sage: K.ring_of_integers().characteristic() # implicit doctest
            0
        """
        if self == 0: return 1
        else: return sage.rings.infinity.infinity

    cdef bint is_rational_c(self):
        return ZZX_deg(self.__numerator) == 0

    def trace(self, K=None):
        """
        Return the absolute or relative trace of this number field
        element.

        If K is given then K must be a subfield of the parent L of self, in
        which case the trace is the relative trace from L to K. In all
        other cases, the trace is the absolute trace down to QQ.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 -132/7*x^2 + x + 1); K
            Number Field in a with defining polynomial x^3 - 132/7*x^2 + x + 1
            sage: a.trace()
            132/7
            sage: (a+1).trace() == a.trace() + 3
            True

        If we are in an order, the trace is an integer::

            sage: K.<zeta> = CyclotomicField(17)
            sage: R = K.ring_of_integers()
            sage: R(zeta).trace().parent()
            Integer Ring

        TESTS::

            sage: F.<z> = CyclotomicField(5) ; t = 3*z**3 + 4*z**2 + 2
            sage: t.trace(F)
            3*z^3 + 4*z^2 + 2
        """
        if K is None:
            trace = self._pari_('x').trace()
            return QQ(trace) if self._parent.is_field() else ZZ(trace)
        return self.matrix(K).trace()

    def norm(self, K=None):
        """
        Return the absolute or relative norm of this number field element.

        If K is given then K must be a subfield of the parent L of self, in
        which case the norm is the relative norm from L to K. In all other
        cases, the norm is the absolute norm down to QQ.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + x^2 + x - 132/7); K
            Number Field in a with defining polynomial x^3 + x^2 + x - 132/7
            sage: a.norm()
            132/7
            sage: factor(a.norm())
            2^2 * 3 * 7^-1 * 11
            sage: K(0).norm()
            0

        Some complicated relatives norms in a tower of number fields.

        ::

            sage: K.<a,b,c> = NumberField([x^2 + 1, x^2 + 3, x^2 + 5])
            sage: L = K.base_field(); M = L.base_field()
            sage: a.norm()
            1
            sage: a.norm(L)
            1
            sage: a.norm(M)
            1
            sage: a
            a
            sage: (a+b+c).norm()
            121
            sage: (a+b+c).norm(L)
            2*c*b - 7
            sage: (a+b+c).norm(M)
            -11

        We illustrate that norm is compatible with towers::

            sage: z = (a+b+c).norm(L); z.norm(M)
            -11

        If we are in an order, the norm is an integer::

            sage: K.<a> = NumberField(x^3-2)
            sage: a.norm().parent()
            Rational Field
            sage: R = K.ring_of_integers()
            sage: R(a).norm().parent()
            Integer Ring

        TESTS::

            sage: F.<z> = CyclotomicField(5)
            sage: t = 3*z**3 + 4*z**2 + 2
            sage: t.norm(F)
            3*z^3 + 4*z^2 + 2
        """
        if K is None:
            norm = self._pari_('x').norm()
            return QQ(norm) if self._parent.is_field() else ZZ(norm)
        return self.matrix(K).determinant()

    def vector(self):
        """
        Return vector representation of self in terms of the basis for the
        ambient number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: (2/3*a - 5/6).vector()
            (-5/6, 2/3)
            sage: (-5/6, 2/3)
            (-5/6, 2/3)
            sage: O = K.order(2*a)
            sage: (O.1).vector()
            (0, 2)
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: (a + b).vector()
            (b, 1)
            sage: O = K.order([a,b])
            sage: (O.1).vector()
            (-b, 1)
            sage: (O.2).vector()
            (1, -b)
        """
        return self.number_field().relative_vector_space()[2](self)

    def charpoly(self, var='x'):
        raise NotImplementedError, "Subclasses of NumberFieldElement must override charpoly()"

    def minpoly(self, var='x'):
        """
        Return the minimal polynomial of this number field element.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+3)
            sage: a.minpoly('x')
            x^2 + 3
            sage: R.<X> = K['X']
            sage: L.<b> = K.extension(X^2-(22 + a))
            sage: b.minpoly('t')
            t^2 - a - 22
            sage: b.absolute_minpoly('t')
            t^4 - 44*t^2 + 487
            sage: b^2 - (22+a)
            0
        """
        return self.charpoly(var).radical() # square free part of charpoly

    def is_integral(self):
        r"""
        Determine if a number is in the ring of integers of this number
        field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: a.is_integral()
            True
            sage: t = (1+a)/2
            sage: t.is_integral()
            True
            sage: t.minpoly()
            x^2 - x + 6
            sage: t = a/2
            sage: t.is_integral()
            False
            sage: t.minpoly()
            x^2 + 23/4

        An example in a relative extension::

            sage: K.<a,b> = NumberField([x^2+1, x^2+3])
            sage: (a+b).is_integral()
            True
            sage: ((a-b)/2).is_integral()
            False
        """
        return all([a in ZZ for a in self.absolute_minpoly()])

    def matrix(self, base=None):
        r"""
        If base is None, return the matrix of right multiplication by the
        element on the power basis `1, x, x^2, \ldots, x^{d-1}` for
        the number field. Thus the *rows* of this matrix give the images of
        each of the `x^i`.

        If base is not None, then base must be either a field that embeds
        in the parent of self or a morphism to the parent of self, in which
        case this function returns the matrix of multiplication by self on
        the power basis, where we view the parent field as a field over
        base.

        INPUT:


        -  ``base`` - field or morphism


        EXAMPLES:

        Regular number field::

            sage: K.<a> = NumberField(QQ['x'].0^3 - 5)
            sage: M = a.matrix(); M
            [0 1 0]
            [0 0 1]
            [5 0 0]
            sage: M.base_ring() is QQ
            True

        Relative number field::

            sage: L.<b> = K.extension(K['x'].0^2 - 2)
            sage: M = b.matrix(); M
            [0 1]
            [2 0]
            sage: M.base_ring() is K
            True

        Absolute number field::

            sage: M = L.absolute_field('c').gen().matrix(); M
            [  0   1   0   0   0   0]
            [  0   0   1   0   0   0]
            [  0   0   0   1   0   0]
            [  0   0   0   0   1   0]
            [  0   0   0   0   0   1]
            [-17 -60 -12 -10   6   0]
            sage: M.base_ring() is QQ
            True

        More complicated relative number field::

            sage: L.<b> = K.extension(K['x'].0^2 - a); L
            Number Field in b with defining polynomial x^2 - a over its base field
            sage: M = b.matrix(); M
            [0 1]
            [a 0]
            sage: M.base_ring() is K
            True

        An example where we explicitly give the subfield or the embedding::

            sage: K.<a> = NumberField(x^4 + 1); L.<a2> = NumberField(x^2 + 1)
            sage: a.matrix(L)
            [ 0  1]
            [a2  0]

        Notice that if we compute all embeddings and choose a different
        one, then the matrix is changed as it should be::

            sage: v = L.embeddings(K)
            sage: a.matrix(v[1])
            [  0   1]
            [-a2   0]

        The norm is also changed::

            sage: a.norm(v[1])
            a2
            sage: a.norm(v[0])
            -a2

        TESTS::

            sage: F.<z> = CyclotomicField(5) ; t = 3*z**3 + 4*z**2 + 2
            sage: t.matrix(F)
            [3*z^3 + 4*z^2 + 2]
        """
        if base is not None:
            if number_field.is_NumberField(base):
                return self._matrix_over_base(base)
            else:
                return self._matrix_over_base_morphism(base)
        # Multiply each power of field generator on
        # the left by this element; make matrix
        # whose rows are the coefficients of the result,
        # and transpose.
        if self.__matrix is None:
            K = self.number_field()
            v = []
            x = K.gen()
            a = K(1)
            d = K.relative_degree()
            for n in range(d):
                v += (a*self).list()
                a *= x
            k = K.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(k, d)
            self.__matrix = M(v)
        return self.__matrix

    def valuation(self, P):
        """
        Returns the valuation of self at a given prime ideal P.

        INPUT:


        -  ``P`` - a prime ideal of the parent of self


        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: b = a^2 + 30
            sage: b.valuation(P)
            1
            sage: type(b.valuation(P))
            <type 'sage.rings.integer.Integer'>
        """
        from number_field_ideal import is_NumberFieldIdeal
        from sage.rings.infinity import infinity
        if not is_NumberFieldIdeal(P):
            if is_NumberFieldElement(P):
                P = self.number_field().fractional_ideal(P)
            else:
                raise TypeError, "P must be an ideal"
        if not P.is_prime():
            # We always check this because it caches the pari prime representation of this ideal.
            raise ValueError, "P must be prime"
        if self == 0:
            return infinity
        return Integer_sage(self.number_field()._pari_().elementval(self._pari_(), P._pari_prime))

    def support(self):
        """
        Return the support of this number field element.

        OUTPUT: A sorted list of the primes ideals at which this number
        field element has nonzero valuation. An error is raised if the
        element is zero.

        EXAMPLES::

            sage: x = ZZ['x'].gen()
            sage: F.<t> = NumberField(x^3 - 2)

        ::

            sage: P5s = F(5).support()
            sage: P5s
            [Fractional ideal (-t^2 - 1), Fractional ideal (t^2 - 2*t - 1)]
            sage: all(5 in P5 for P5 in P5s)
            True
            sage: all(P5.is_prime() for P5 in P5s)
            True
            sage: [ P5.norm() for P5 in P5s ]
            [5, 25]

        TESTS:

        It doesn't make sense to factor the ideal (0)::

            sage: F(0).support()
            Traceback (most recent call last):
            ...
            ArithmeticError: Support of 0 is not defined.
        """
        if self.is_zero():
            raise ArithmeticError, "Support of 0 is not defined."
        return self.number_field().primes_above(self)

    def _matrix_over_base(self, L):
        """
        Return the matrix of self over the base field L.

        EXAMPLES::

            sage: K.<a> = NumberField(ZZ['x'].0^3-2, 'a')
            sage: L.<b> = K.extension(ZZ['x'].0^2+3, 'b')
            sage: L(a)._matrix_over_base(K) == L(a).matrix()
            True
        """
        K = self.number_field()
        E = L.embeddings(K)
        if len(E) == 0:
            raise ValueError, "no way to embed L into parent's base ring K"
        phi = E[0]
        return self._matrix_over_base_morphism(phi)

    def _matrix_over_base_morphism(self, phi):
        """
        Return the matrix of self over a specified base, where phi gives a
        map from the specified base to self.parent().

        EXAMPLES::

            sage: F.<alpha> = NumberField(ZZ['x'].0^5-2)
            sage: h = Hom(QQ,F)([1])
            sage: alpha._matrix_over_base_morphism(h) == alpha.matrix()
            True
            sage: alpha._matrix_over_base_morphism(h) == alpha.matrix(QQ)
            True
        """
        L = phi.domain()

        ## the code below doesn't work if the morphism is
        ## over QQ, since QQ.primitive_element() doesn't
        ## make sense
        if L is QQ:
            K = phi.codomain()
            if K != self.number_field():
                raise ValueError, "codomain of phi must be parent of self"
            ## the variable name is irrelevant below, because the
            ## matrix is over QQ
            F = K.absolute_field('alpha')
            from_f, to_F = F.structure()
            return to_F(self).matrix()

        alpha = L.primitive_element()
        beta = phi(alpha)
        K = phi.codomain()
        if K != self.number_field():
            raise ValueError, "codomain of phi must be parent of self"

        # Construct a relative extension over L (= QQ(beta))
        M = K.relativize(beta, ('a','b'))
                     # variable name a is OK, since this is temporary

        # Carry self over to M.
        from_M, to_M = M.structure()
        try:
            z = to_M(self)
        except Exception:
            return to_M, self, K, beta

        # Compute the relative matrix of self, but in M
        R = z.matrix()

        # Map back to L.
        psi = M.base_field().hom([alpha])
        return R.apply_morphism(psi)


    def list(self):
        """
        Return the list of coefficients of self written in terms of a power
        basis.
        """
        # Power basis list is total nonsense, unless the parent of self is an
        # absolute extension.
        raise NotImplementedError


cdef class NumberFieldElement_absolute(NumberFieldElement):

    def _pari_(self, var='x'):
        """
        Return PARI C-library object corresponding to self.

        EXAMPLES::

            sage: k.<j> = QuadraticField(-1)
            sage: j._pari_('j')
            Mod(j, j^2 + 1)
            sage: pari(j)
            Mod(x, x^2 + 1)

        ::

            sage: y = QQ['y'].gen()
            sage: k.<j> = NumberField(y^3 - 2)
            sage: pari(j)
            Mod(x, x^3 - 2)

        By default the variable name is 'x', since in PARI many variable
        names are reserved::

            sage: theta = polygen(QQ, 'theta')
            sage: M.<theta> = NumberField(theta^2 + 1)
            sage: pari(theta)
            Mod(x, x^2 + 1)

        If you try do coerce a generator called I to PARI, hell may break
        loose::

            sage: k.<I> = QuadraticField(-1)
            sage: I._pari_('I')
            Traceback (most recent call last):
            ...
            PariError: forbidden (45)

        Instead, request the variable be named different for the coercion::

            sage: pari(I)
            Mod(x, x^2 + 1)
            sage: I._pari_('i')
            Mod(i, i^2 + 1)
            sage: I._pari_('II')
            Mod(II, II^2 + 1)
        """
        try:
            return self.__pari[var]
        except KeyError:
            pass
        except TypeError:
            self.__pari = {}
        if var is None:
            var = self.number_field().variable_name()
        f = self.polynomial()._pari_()
        gp = self.number_field().polynomial()
        if gp.variable_name() != 'x':
            gp = gp.change_variable_name('x')
        g = gp._pari_()
        gv = str(gp.parent().gen())
        if var != 'x':
            f = f.subst("x",var)
        if var != gv:
            g = g.subst(gv, var)
        h = f.Mod(g)
        self.__pari[var] = h
        return h

    def _magma_init_(self, magma):
        """
        Return Magma version of this number field element.

        INPUT:


        -  ``magma`` - a Magma interpreter


        OUTPUT: MagmaElement that has parent the Magma object corresponding
        to the parent number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: a._magma_init_(magma)            # optional - magma
            '(_sage_[...]![0, 1, 0])'
            sage: magma((2/3)*a^2 - 17/3)          # optional - magma
            1/3*(2*a^2 - 17)

        An element of a cyclotomic field.

        ::

            sage: K = CyclotomicField(9)
            sage: K.gen()
            zeta9
            sage: K.gen()._magma_init_(magma)     # optional - magma
            '(_sage_[...]![0, 1, 0, 0, 0, 0])'
            sage: magma(K.gen())                  # optional - magma
            zeta9
        """
        K = magma(self.parent())
        return '(%s!%s)'%(K.name(), self.list())

    cdef void _parent_poly_c_(self, ZZX_c *num, ZZ_c *den):
        """
        I believe this function should be removed since I've put the
        pointer __fld_numerator and __fld_denominator in the element
        class. I'm not going to remove it quite yet, but feel free to
        remove it if you agree with me that it should go.
        """
        cdef ntl_ZZX _num
        cdef ntl_ZZ _den
        _num, _den = self.number_field().polynomial_ntl()
        num[0] = _num.x
        den[0] = _den.x

    def absolute_charpoly(self, var='x', algorithm=None):
        r"""
        Return the characteristic polynomial of this element over

        For the meaning of the optional argument algorithm, see :meth:`charpoly`.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a> = NumberField(x^4 + 2, 'a')
            sage: a.absolute_charpoly()
            x^4 + 2
            sage: a.absolute_charpoly('y')
            y^4 + 2
            sage: (-a^2).absolute_charpoly()
            x^4 + 4*x^2 + 4
            sage: (-a^2).absolute_minpoly()
            x^2 + 2

            sage: a.absolute_charpoly(algorithm='pari') == a.absolute_charpoly(algorithm='sage')
            True
        """
        # this hack is necessary because quadratic fields override
        # charpoly(), and they don't take the argument 'algorithm'
        if algorithm is None:
            return self.charpoly(var)
        return self.charpoly(var, algorithm)

    def absolute_minpoly(self, var='x', algorithm=None):
        r"""
        Return the minimal polynomial of this element over
        `\QQ`.

        For the meaning of the optional argument algorithm, see :meth:`charpoly`.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: f = x^10 - 5*x^9 + 15*x^8 - 68*x^7 + 81*x^6 - 221*x^5 + 141*x^4 - 242*x^3 - 13*x^2 - 33*x - 135
            sage: K.<a> = NumberField(f, 'a')
            sage: a.absolute_charpoly()
            x^10 - 5*x^9 + 15*x^8 - 68*x^7 + 81*x^6 - 221*x^5 + 141*x^4 - 242*x^3 - 13*x^2 - 33*x - 135
            sage: a.absolute_charpoly('y')
            y^10 - 5*y^9 + 15*y^8 - 68*y^7 + 81*y^6 - 221*y^5 + 141*y^4 - 242*y^3 - 13*y^2 - 33*y - 135
            sage: b = -79/9995*a^9 + 52/9995*a^8 + 271/9995*a^7 + 1663/9995*a^6 + 13204/9995*a^5 + 5573/9995*a^4 + 8435/1999*a^3 - 3116/9995*a^2 + 7734/1999*a + 1620/1999
            sage: b.absolute_charpoly()
            x^10 + 10*x^9 + 25*x^8 - 80*x^7 - 438*x^6 + 80*x^5 + 2950*x^4 + 1520*x^3 - 10439*x^2 - 5130*x + 18225
            sage: b.absolute_minpoly()
            x^5 + 5*x^4 - 40*x^2 - 19*x + 135

            sage: b.absolute_minpoly(algorithm='pari') == b.absolute_minpoly(algorithm='sage')
            True
        """
        # this hack is necessary because quadratic fields override
        # minpoly(), and they don't take the argument 'algorithm'
        if algorithm is None:
            return self.minpoly(var)
        return self.minpoly(var, algorithm)

    def charpoly(self, var='x', algorithm=None):
        r"""
        The characteristic polynomial of this element, over
        `\QQ` if self is an element of a field, and over
        `\ZZ` is self is an element of an order.

        This is the same as ``self.absolute_charpoly`` since
        this is an element of an absolute extension.

        The optional argument algorithm controls how the
        characteristic polynomial is computed: 'pari' uses Pari,
        'sage' uses charpoly for Sage matrices.  The default value
        None means that 'pari' is used for small degrees (up to the
        value of the constant TUNE_CHARPOLY_NF, currently at 25),
        otherwise 'sage' is used.  The constant TUNE_CHARPOLY_NF
        should give reasonable performance on all architectures;
        however, if you feel the need to customize it to your own
        machine, see trac ticket 5213 for a tuning script.

        EXAMPLES:

        We compute the characteristic polynomial of the cube root of `2`.

        ::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: a.charpoly('x')
            x^3 - 2
            sage: a.charpoly('y').parent()
            Univariate Polynomial Ring in y over Rational Field

        TESTS::

            sage: R = K.ring_of_integers()
            sage: R(a).charpoly()
            x^3 - 2
            sage: R(a).charpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring

            sage: R(a).charpoly(algorithm='pari') == R(a).charpoly(algorithm='sage')
            True
        """
        if algorithm is None:
            if self._parent.degree() <= TUNE_CHARPOLY_NF:
                algorithm = 'pari'
            else:
                algorithm = 'sage'
        R = self._parent.base_ring()[var]
        if algorithm == 'pari':
            return R(self._pari_('x').charpoly())
        if algorithm == 'sage':
            return R(self.matrix().charpoly())

    def minpoly(self, var='x', algorithm=None):
        """
        Return the minimal polynomial of this number field element.

        For the meaning of the optional argument algorithm, see charpoly().

        EXAMPLES:

        We compute the characteristic polynomial of cube root of `2`.

        ::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: a.minpoly('x')
            x^3 - 2
            sage: a.minpoly('y').parent()
            Univariate Polynomial Ring in y over Rational Field

        TESTS::

            sage: R = K.ring_of_integers()
            sage: R(a).minpoly()
            x^3 - 2
            sage: R(a).minpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring

            sage: R(a).minpoly(algorithm='pari') == R(a).minpoly(algorithm='sage')
            True

        """
        return self.charpoly(var, algorithm).radical() # square free part of charpoly

    def list(self):
        """
        Return the list of coefficients of self written in terms of a power
        basis.

        EXAMPLE::

            sage: K.<z> = CyclotomicField(3)
            sage: (2+3/5*z).list()
            [2, 3/5]
            sage: (5*z).list()
            [0, 5]
            sage: K(3).list()
            [3, 0]
        """
        n = self.number_field().degree()
        v = self._coefficients()
        z = sage.rings.rational.Rational(0)
        return v + [z]*(n - len(v))

    def inverse_mod(self, I):
        """
        Returns the inverse of self mod the integral ideal I.

        INPUT:

        -  ``I`` - may be an ideal of self.parent(), or an element or list
           of elements of self.parent() generating a nonzero ideal. A TypeError
           is raised if I is non-integral, and a ValueError if the generators
           are all zero. A ZeroDivisionError is raised if I + (x) != (1).

        NOTE: It's not implemented yet for non-integral elements.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3)
            sage: d = 3*a + 1
            sage: d.inverse_mod(N)
            1
        """
        R = I.number_field().ring_of_integers()
        try:
            return I.small_residue(_inverse_mod_generic(R(self), I))
        except TypeError:
            raise NotImplementedError, "inverse_mod is not implemented for non-integral elements"


cdef class NumberFieldElement_relative(NumberFieldElement):
    r"""
    The current relative number field element implementation
    does everything in terms of absolute polynomials.

    All conversions from relative polynomials, lists, vectors, etc
    should happen in the parent.
    """
    def __init__(self, parent, f):
        NumberFieldElement.__init__(self, parent, f)

    def __getitem__(self, n):
        """
        Return the n-th coefficient of this relative number field element, written
        as a polynomial in the generator.

        Note that `n` must be between 0 and `d-1`, where
        `d` is the relative degree of the number field.

        EXAMPLES::

            sage: K.<a, b> = NumberField([x^3 - 5, x^2 + 3])
            sage: c = (a + b)^3; c
            3*b*a^2 - 9*a - 3*b + 5
            sage: c[0]
            -3*b + 5

        We illustrate bounds checking::

            sage: c[-1]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and the relative degree minus 1.
            sage: c[4]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and the relative degree minus 1.

        The list method implicitly calls ``__getitem__``::

            sage: list(c)
            [-3*b + 5, -9, 3*b]
            sage: K(list(c)) == c
            True
        """
        if n < 0 or n >= self.parent().relative_degree():
            raise IndexError, "index must be between 0 and the relative degree minus 1."
        return self.vector()[n]

    def list(self):
        """
        Return the list of coefficients of self written in terms of a power
        basis.

        EXAMPLES::

            sage: K.<a,b> = NumberField([x^3+2, x^2+1])
            sage: a.list()
            [0, 1, 0]
            sage: v = (K.base_field().0 + a)^2 ; v
            a^2 + 2*b*a - 1
            sage: v.list()
            [-1, 2*b, 1]
        """
        return self.vector().list()

    def _pari_(self, var='x'):
        """
        Return PARI C-library object corresponding to self.

        EXAMPLES:

        By default the variable name is 'x', since in PARI many
        variable names are reserved.

        ::

            sage: y = QQ['y'].gen()
            sage: k.<j> = NumberField([y^2 - 7, y^3 - 2])
            sage: pari(j)
            Mod(42/5515*x^5 - 9/11030*x^4 - 196/1103*x^3 + 273/5515*x^2 + 10281/5515*x + 4459/11030, x^6 - 21*x^4 + 4*x^3 + 147*x^2 + 84*x - 339)
            sage: j^2
            7
            sage: pari(j)^2
            Mod(7, x^6 - 21*x^4 + 4*x^3 + 147*x^2 + 84*x - 339)
            sage: (j^2)._pari_('y')
            Mod(7, y^6 - 21*y^4 + 4*y^3 + 147*y^2 + 84*y - 339)

            sage: K.<a> = NumberField(x^2 + 2, 'a')
            sage: K(1)._pari_()
            Mod(1, x^2 + 2)
            sage: K(1)._pari_('t')
            Mod(1, t^2 + 2)

            sage: K.gen()._pari_()
            Mod(x, x^2 + 2)
            sage: K.gen()._pari_('t')
            Mod(t, t^2 + 2)

            At this time all elements, even relative elements, are
            represented as absolute polynomials:

            sage: K.<a> = NumberField(x^2 + 2, 'a')
            sage: L.<b> = NumberField(K['x'].0^2 + a, 'b')
            sage: L(1)._pari_()
            Mod(1, x^4 + 2)
            sage: L(1)._pari_('t')
            Mod(1, t^4 + 2)
            sage: L.gen()._pari_()
            Mod(x, x^4 + 2)
            sage: L.gen()._pari_('t')
            Mod(t, t^4 + 2)

            sage: M.<c> = NumberField(L['x'].0^3 + b, 'c')
            sage: M(1)._pari_()
            Mod(1, x^12 + 2)
            sage: M(1)._pari_('t')
            Mod(1, t^12 + 2)
            sage: M.gen()._pari_()
            Mod(x, x^12 + 2)
            sage: M.gen()._pari_('t')
            Mod(t, t^12 + 2)
        """
        try:
            return self.__pari[var]
        except KeyError:
            pass
        except TypeError:
            self.__pari = {}
        g = self.parent().pari_polynomial()
        f = self.polynomial()._pari_()
        f = f.subst('x', var)
        g = g.subst('x', var)
        h = f.Mod(g)
        h = f.Mod(g)
        self.__pari[var] = h
        return h

    cdef void _parent_poly_c_(self, ZZX_c *num, ZZ_c *den):
        """
        I believe this function should be removed since I've put the
        pointer __fld_numerator and __fld_denominator in the element
        class. I'm not going to remove it quite yet, but feel free to
        remove it if you agree with me that it should go.
        """
        f = self.number_field().absolute_polynomial()
        _ntl_poly(f, num, den)

    def __repr__(self):
        K = self.number_field()
        # Compute representation of self in terms of relative vector space.
        R = K.base_field()[K.variable_name()]
        return repr(R(self.list()))

    def _latex_(self):
        r"""
        Returns the latex representation for this element.

        EXAMPLES::

            sage: C.<zeta> = CyclotomicField(12)
            sage: PC.<x> = PolynomialRing(C)
            sage: K.<alpha> = NumberField(x^2 - 7)
            sage: latex((alpha + zeta)^4)
            \left(4 \zeta_{12}^{3} + 28 \zeta_{12}\right) \alpha + 43 \zeta_{12}^{2} + 48
            sage: PK.<y> = PolynomialRing(K)
            sage: L.<beta> = NumberField(y^3 + y + alpha)
            sage: latex((beta + zeta)^3)
            3 \zeta_{12} \beta^{2} + \left(3 \zeta_{12}^{2} - 1\right) \beta - \alpha + \zeta_{12}^{3}
        """
        K = self.number_field()
        R = K.base_field()[K.variable_name()]
        return R(self.list())._latex_()

    def charpoly(self, var='x'):
        r"""
        The characteristic polynomial of this element over its base field.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a, b> = QQ.extension([x^2 + 2, x^5 + 400*x^4 + 11*x^2 + 2])
            sage: a.charpoly()
            x^2 + 2
            sage: b.charpoly()
            x^2 - 2*b*x + b^2
            sage: b.minpoly()
            x - b

            sage: K.<a, b> = NumberField([x^2 + 2, x^2 + 1000*x + 1])
            sage: y = K['y'].0
            sage: L.<c> = K.extension(y^2 + a*y + b)
            sage: c.charpoly()
            x^2 + a*x + b
            sage: c.minpoly()
            x^2 + a*x + b
            sage: L(a).charpoly()
            x^2 - 2*a*x - 2
            sage: L(a).minpoly()
            x - a
            sage: L(b).charpoly()
            x^2 - 2*b*x - 1000*b - 1
            sage: L(b).minpoly()
            x - b
        """
        R = self._parent.base_ring()[var]
        return R(self.matrix().charpoly())

    def absolute_charpoly(self, var='x', algorithm=None):
        r"""
        The characteristic polynomial of this element over
        `\QQ`.

        We construct a relative extension and find the characteristic
        polynomial over `\QQ`.

        The optional argument algorithm controls how the
        characteristic polynomial is computed: 'pari' uses Pari,
        'sage' uses charpoly for Sage matrices.  The default value
        None means that 'pari' is used for small degrees (up to the
        value of the constant TUNE_CHARPOLY_NF, currently at 25),
        otherwise 'sage' is used.  The constant TUNE_CHARPOLY_NF
        should give reasonable performance on all architectures;
        however, if you feel the need to customize it to your own
        machine, see trac ticket 5213 for a tuning script.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: S.<X> = K[]
            sage: L.<b> = NumberField(X^3 + 17); L
            Number Field in b with defining polynomial X^3 + 17 over its base field
            sage: b.absolute_charpoly()
            x^9 + 51*x^6 + 867*x^3 + 4913
            sage: b.charpoly()(b)
            0
            sage: a = L.0; a
            b
            sage: a.absolute_charpoly('x')
            x^9 + 51*x^6 + 867*x^3 + 4913
            sage: a.absolute_charpoly('y')
            y^9 + 51*y^6 + 867*y^3 + 4913

            sage: a.absolute_charpoly(algorithm='pari') == a.absolute_charpoly(algorithm='sage')
            True
        """
        if algorithm is None:
            # this might not be the optimal condition; maybe it should
            # be .degree() instead of .absolute_degree()
            # there are too many bugs in relative number fields to
            # figure this out now
            if self._parent.absolute_degree() <= TUNE_CHARPOLY_NF:
                algorithm = 'pari'
            else:
                algorithm = 'sage'
        g = self.polynomial()  # in QQ[x]
        R = QQ[var]
        if algorithm == 'pari':
            f = self.number_field().pari_polynomial()  # # field is QQ[x]/(f)
            return R((g._pari_().Mod(f)).charpoly()).change_variable_name(var)
        if algorithm == 'sage':
            return R(self.matrix(QQ).charpoly())

    def absolute_minpoly(self, var='x', algorithm=None):
        r"""
        Return the minimal polynomial over `\QQ` of this element.

        For the meaning of the optional argument algorithm, see :meth:`absolute_charpoly`.

        EXAMPLES::

            sage: K.<a, b> = NumberField([x^2 + 2, x^2 + 1000*x + 1])
            sage: y = K['y'].0
            sage: L.<c> = K.extension(y^2 + a*y + b)
            sage: c.absolute_charpoly()
            x^8 - 1996*x^6 + 996006*x^4 + 1997996*x^2 + 1
            sage: c.absolute_minpoly()
            x^8 - 1996*x^6 + 996006*x^4 + 1997996*x^2 + 1
            sage: L(a).absolute_charpoly()
            x^8 + 8*x^6 + 24*x^4 + 32*x^2 + 16
            sage: L(a).absolute_minpoly()
            x^2 + 2
            sage: L(b).absolute_charpoly()
            x^8 + 4000*x^7 + 6000004*x^6 + 4000012000*x^5 + 1000012000006*x^4 + 4000012000*x^3 + 6000004*x^2 + 4000*x + 1
            sage: L(b).absolute_minpoly()
            x^2 + 1000*x + 1
        """
        return self.absolute_charpoly(var, algorithm).radical()

    def valuation(self, P):
        """
        Returns the valuation of self at a given prime ideal P.

        INPUT:


        -  ``P`` - a prime ideal of relative number field which is the parent of self


        EXAMPLES::

            sage: K.<a, b, c> = NumberField([x^2 - 2, x^2 - 3, x^2 - 5])
            sage: P = K.prime_factors(5)[0]
            sage: (2*a + b - c).valuation(P)
            1
        """
        P_abs = P.absolute_ideal()
        abs = P_abs.number_field()
        to_abs = abs.structure()[1]
        return to_abs(self).valuation(P_abs)


cdef class OrderElement_absolute(NumberFieldElement_absolute):
    """
    Element of an order in an absolute number field.

    EXAMPLES::

        sage: K.<a> = NumberField(x^2 + 1)
        sage: O2 = K.order(2*a)
        sage: w = O2.1; w
        2*a
        sage: parent(w)
        Order in Number Field in a with defining polynomial x^2 + 1

        sage: w.absolute_charpoly()
        x^2 + 4
        sage: w.absolute_charpoly().parent()
        Univariate Polynomial Ring in x over Integer Ring
        sage: w.absolute_minpoly()
        x^2 + 4
        sage: w.absolute_minpoly().parent()
        Univariate Polynomial Ring in x over Integer Ring
    """
    def __init__(self, order, f):
        K = order.number_field()
        NumberFieldElement_absolute.__init__(self, K, f)
        self._number_field = K
        (<Element>self)._parent = order

    cdef _new(self):
        """
        Quickly creates a new initialized NumberFieldElement with the same
        parent as self.

        EXAMPLES:

        This is called implicitly in multiplication::

            sage: O = EquationOrder(x^3 + 18, 'a')
            sage: O.1 * O.1 * O.1
            -18
        """
        cdef OrderElement_absolute x
        x = <OrderElement_absolute>PY_NEW_SAME_TYPE(self)
        x._parent = self._parent
        x._number_field = self._parent.number_field()
        x.__fld_numerator = self.__fld_numerator
        x.__fld_denominator = self.__fld_denominator
        return x

    cdef number_field(self):
        return self._number_field

    cpdef RingElement _div_(self, RingElement other):
        r"""
        Implement division, checking that the result has the right parent.
        It's not so crucial what the parent actually is, but it is crucial
        that the returned value really is an element of its supposed
        parent! This fixes trac #4190.

        EXAMPLES::

            sage: K = NumberField(x^3 - 17, 'a')
            sage: OK = K.ring_of_integers()
            sage: a = OK(K.gen())
            sage: (17/a) in OK
            True
            sage: (17/a).parent() is K
            True
            sage: (17/(2*a)).parent() is K
            True
            sage: (17/(2*a)) in OK
            False
        """
        cdef NumberFieldElement_absolute x
        x = NumberFieldElement_absolute._div_(self, other)
        return self._parent.number_field()(x)

    def inverse_mod(self, I):
        r"""
        Return an inverse of self modulo the given ideal.

        INPUT:


        -  ``I`` - may be an ideal of self.parent(), or an
           element or list of elements of self.parent() generating a nonzero
           ideal. A TypeError is raised if I is non-integral, and a ValueError
           if the generators are all zero. A ZeroDivisionError is raised if I
           + (x) != (1).


        EXAMPLES::

            sage: OE = NumberField(x^3 - x + 2, 'w').ring_of_integers()
            sage: w = OE.ring_generators()[0]
            sage: w.inverse_mod(13*OE)
            -7*w^2 - 13*w + 7
            sage: w * (w.inverse_mod(13)) - 1 in 13*OE
            True
            sage: w.inverse_mod(2*OE)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: w is not invertible modulo Fractional ideal (2)
        """
        return _inverse_mod_generic(self, I)

    def __invert__(self):
        r"""
        Implement inversion, checking that the return value has the right
        parent. See trac #4190.

        EXAMPLE::

            sage: K = NumberField(x^3 -x + 2, 'a')
            sage: OK = K.ring_of_integers()
            sage: a = OK(K.gen())
            sage: (~a).parent() is K
            True
            sage: (~a) in OK
            False
            sage: a**(-1) in OK
            False
        """
        return self._parent.number_field()(NumberFieldElement_absolute.__invert__(self))

cdef class OrderElement_relative(NumberFieldElement_relative):
    """
    Element of an order in a relative number field.

    EXAMPLES::

        sage: O = EquationOrder([x^2 + x + 1, x^3 - 2],'a,b')
        sage: c = O.1; c
        (-2*b^2 - 2)*a - 2*b^2 - b
        sage: type(c)
        <type 'sage.rings.number_field.number_field_element.OrderElement_relative'>
    """
    def __init__(self, order, f):
        K = order.number_field()
        NumberFieldElement_relative.__init__(self, K, f)
        (<Element>self)._parent = order
        self._number_field = K

    cdef number_field(self):
        return self._number_field

    cdef _new(self):
        """
        Quickly creates a new initialized NumberFieldElement with the same
        parent as self.

        EXAMPLES:

        This is called implicitly in multiplication::

            sage: O = EquationOrder([x^2 + 18, x^3 + 2], 'a,b')
            sage: c = O.1 * O.2; c
            (-23321*b^2 - 9504*b + 10830)*a + 10152*b^2 - 104562*b - 110158
            sage: parent(c) == O
            True
        """
        cdef OrderElement_relative x
        x = <OrderElement_relative>PY_NEW_SAME_TYPE(self)
        x._parent = self._parent
        x._number_field = self._parent.number_field()
        x.__fld_numerator = self.__fld_numerator
        x.__fld_denominator = self.__fld_denominator
        return x

    cpdef RingElement _div_(self, RingElement other):
        r"""
        Implement division, checking that the result has the right parent.
        It's not so crucial what the parent actually is, but it is crucial
        that the returned value really is an element of its supposed
        parent. This fixes trac #4190.

        EXAMPLES::

            sage: K1.<a> = NumberField(x^3 - 17)
            sage: R.<y> = K1[]
            sage: K2 = K1.extension(y^2 - a, 'b')
            sage: OK2 = K2.order(K2.gen()) # (not maximal)
            sage: b = OK2.gens()[1]; b
            b
            sage: (17/b).parent() is K2
            True
            sage: (17/b) in OK2 # not implemented (#4193)
            True
            sage: (17/b^7) in OK2
            False
        """
        cdef NumberFieldElement_relative x
        x = NumberFieldElement_relative._div_(self, other)
        return self._parent.number_field()(x)

    def __invert__(self):
        r"""
        Implement division, checking that the result has the right parent.
        See trac #4190.

        EXAMPLES::

            sage: K1.<a> = NumberField(x^3 - 17)
            sage: R.<y> = K1[]
            sage: K2 = K1.extension(y^2 - a, 'b')
            sage: OK2 = K2.order(K2.gen()) # (not maximal)
            sage: b = OK2.gens()[1]; b
            b
            sage: b.parent() is OK2
            True
            sage: (~b).parent() is K2
            True
            sage: (~b) in OK2 # not implemented (#4193)
            False
            sage: b**(-1) in OK2 # not implemented (#4193)
            False
        """
        return self._parent.number_field()(NumberFieldElement_relative.__invert__(self))

    def charpoly(self, var='x'):
        r"""
        The characteristic polynomial of this order element over its base ring.

        This special implementation works around bug \#4738.  At this
        time the base ring of relative order elements is ZZ; it should
        be the ring of integers of the base field.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order(); OK.basis()
            [1, 1/2*a - 1/2*b, -1/2*b*a + 1/2, a]
            sage: charpoly(OK.1)
            x^2 + b*x + 1
            sage: charpoly(OK.1).parent()
            Univariate Polynomial Ring in x over Maximal Order in Number Field in b with defining polynomial x^2 - 3
            sage: [ charpoly(t) for t in OK.basis() ]
            [x^2 - 2*x + 1, x^2 + b*x + 1, x^2 - x + 1, x^2 + 1]
        """
        R = self.parent().number_field().base_field().ring_of_integers()[var]
        return R(self.matrix().charpoly(var))

    def minpoly(self, var='x'):
        r"""
        The minimal polynomial of this order element over its base ring.

        This special implementation works around bug \#4738.  At this
        time the base ring of relative order elements is ZZ; it should
        be the ring of integers of the base field.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order(); OK.basis()
            [1, 1/2*a - 1/2*b, -1/2*b*a + 1/2, a]
            sage: minpoly(OK.1)
            x^2 + b*x + 1
            sage: charpoly(OK.1).parent()
            Univariate Polynomial Ring in x over Maximal Order in Number Field in b with defining polynomial x^2 - 3
            sage: _, u, _, v = OK.basis()
            sage: t = 2*u - v; t
            -b
            sage: t.charpoly()
            x^2 + 2*b*x + 3
            sage: t.minpoly()
            x + b

            sage: t.absolute_charpoly()
            x^4 - 6*x^2 + 9
            sage: t.absolute_minpoly()
            x^2 - 3
        """
        K = self.parent().number_field()
        R = K.base_field().ring_of_integers()[var]
        return R(K(self).minpoly(var))

    def absolute_charpoly(self, var='x'):
        r"""
        The absolute characteristic polynomial of this order element over ZZ.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order()
            sage: _, u, _, v = OK.basis()
            sage: t = 2*u - v; t
            -b
            sage: t.absolute_charpoly()
            x^4 - 6*x^2 + 9
            sage: t.absolute_minpoly()
            x^2 - 3
            sage: t.absolute_charpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring
        """
        K = self.parent().number_field()
        R = ZZ[var]
        return R(K(self).absolute_charpoly(var))

    def absolute_minpoly(self, var='x'):
        r"""
        The absolute minimal polynomial of this order element over ZZ.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order()
            sage: _, u, _, v = OK.basis()
            sage: t = 2*u - v; t
            -b
            sage: t.absolute_charpoly()
            x^4 - 6*x^2 + 9
            sage: t.absolute_minpoly()
            x^2 - 3
            sage: t.absolute_minpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring
        """
        K = self.parent().number_field()
        R = ZZ[var]
        return R(K(self).absolute_minpoly(var))



class CoordinateFunction:
    def __init__(self, NumberFieldElement alpha, W, to_V):
        self.__alpha = alpha
        self.__W = W
        self.__to_V = to_V
        self.__K = alpha.number_field()

    def __repr__(self):
        return "Coordinate function that writes elements in terms of the powers of %s"%self.__alpha

    def alpha(self):
        return self.__alpha

    def __call__(self, x):
        return self.__W.coordinates(self.__to_V(self.__K(x)))




#################

cdef void _ntl_poly(f, ZZX_c *num, ZZ_c *den):
    cdef long i
    cdef ZZ_c coeff
    cdef ntl_ZZX _num
    cdef ntl_ZZ _den

    __den = f.denominator()
    (<Integer>ZZ(__den))._to_ZZ(den)

    __num = f * __den
    for i from 0 <= i <= __num.degree():
        (<Integer>ZZ(__num[i]))._to_ZZ(&coeff)
        ZZX_SetCoeff( num[0], i, coeff )


