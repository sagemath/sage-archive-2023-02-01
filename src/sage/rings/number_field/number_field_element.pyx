"""
Number Field Elements

AUTHORS:
    -- Joel B. Mohler (2007-03-09) - First reimplementation into pyrex
"""

# TODO -- relative extensions need to be completely rewritten, so one
# can get easy access to representation of elements in their relative
# form.  Functions like matrix below can't be done until relative
# extensions are re-written this way.  Also there needs to be class
# hierarchy for number field elements, integers, etc.  This is a
# nontrivial project, and it needs somebody to attack it.  I'm amazed
# how long this has gone unattacked.

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

include "../../ext/stdsage.pxi"

import sage.rings.field_element
import sage.rings.infinity
import sage.rings.polynomial_element
import sage.rings.polynomial_ring
import sage.rings.rational_field
import sage.rings.rational
import sage.rings.integer_ring
import sage.rings.integer
import sage.rings.arith

import sage.rings.number_field.number_field

from sage.libs.ntl.ntl cimport ntl_ZZ, ntl_ZZX
from sage.rings.integer_ring cimport IntegerRing_class

from sage.libs.all import pari_gen
from sage.libs.pari.gen import PariError

QQ = sage.rings.rational_field.QQ
ZZ = sage.rings.integer_ring.ZZ
Integer_sage = sage.rings.integer.Integer

def is_NumberFieldElement(x):
    return isinstance(x, NumberFieldElement)

def __create__NumberFieldElement_version0(parent, poly):
    return NumberFieldElement(parent, poly)

cdef class NumberFieldElement(FieldElement):
    """
    An element of a number field.
    """

    cdef NumberFieldElement _new(self):
        """
        Quickly creates a new initialized NumberFieldElement with the same parent as self.
        """
        cdef NumberFieldElement x
        x = PY_NEW(NumberFieldElement)
        x._parent = self._parent
        return x

    def __init__(self, parent, f):
        """
        INPUT:
            parent -- a number field
            f -- defines an element of a number field.

        EXAMPLES:
        The following examples illustrate creation of elements of
        number fields, and some basic arithmetic.

        First we define a polynomial over Q.
            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^2 + 1

        Next we use f to define the number field.
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
            sage: K.<b> = NumberField(x^3 - 2)
            sage: b = K.gen()
            sage: b^3
            2
            sage: (b^2 + b + 1)^3
            12*b^2 + 15*b + 19

        This example illustrates save and load:
            sage: K.<a> = NumberField(x^17 - 2)
            sage: s = a^15 - 19*a + 3
            sage: loads(s.dumps()) == s
            True
        """
        sage.rings.field_element.FieldElement.__init__(self, parent)

        cdef ntl_c_ZZ coeff
        if isinstance(f, (int, long, Integer_sage)):
            # set it up and exit immediately
            # fast pathway
            (<Integer>ZZ(f))._to_ZZ(&coeff)
            SetCoeff( self.__numerator, 0, coeff )
            conv_ZZ_int( self.__denominator, 1 )
            return

        ppr = parent.polynomial_ring()
        if isinstance(parent, sage.rings.number_field.number_field.NumberField_extension):
            ppr = parent.base_field().polynomial_ring()

        if isinstance(f, pari_gen):
            f = f.lift()
            f = ppr(f)
        if not isinstance(f, sage.rings.polynomial_element.Polynomial):
            f = ppr(f)
        if f.degree() >= parent.degree():
            if isinstance(parent, sage.rings.number_field.number_field.NumberField_extension):
                f %= parent.absolute_polynomial()
            else:
                f %= parent.polynomial()

        # Set Denominator
        den = f.denominator()
        (<Integer>ZZ(den))._to_ZZ(&self.__denominator)

        cdef long i
        num = f * den
        for i from 0 <= i <= num.degree():
            (<Integer>ZZ(num[i]))._to_ZZ(&coeff)
            SetCoeff( self.__numerator, i, coeff )

    def _lift_cyclotomic_element(self, new_parent):
        """
            Creates an element of the passed field from this field.  This is specific to creating elements in a
        cyclotomic field from elements in another cyclotomic field.  This function aims to make this common
        coercion extremely fast!

        EXAMPLES:
            sage: C.<zeta5>=CyclotomicField(5)
            sage: CyclotomicField(10)(zeta5+1)  # The function _lift_cyclotomic_element does the heavy lifting in the background
            zeta10^2 + 1
            sage: (zeta5+1)._lift_cyclotomic_element(CyclotomicField(10))  # There is rarely a purpose to call this function directly
            zeta10^2 + 1

        AUTHOR:
            Joel B. Mohler
        """
        # Right now, I'm a little confused why quadratic extension fields have a zeta_order function
        # I would rather they not have this function since I don't want to do this isinstance check here.
        if not isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_cyclotomic) or not isinstance(new_parent, sage.rings.number_field.number_field.NumberField_cyclotomic):
            raise TypeError, "The field and the new parent field must both be cyclotomic fields."

        try:
            small_order = self.parent().zeta_order()
            large_order = new_parent.zeta_order()
        except AttributeError:
            raise TypeError, "The field and the new parent field must both be cyclotomic fields."

        try:
            _rel = ZZ(large_order / small_order)
        except TypeError:
            raise TypeError, "The zeta_order of the new field must be a multiple of the zeta_order of the original."

        cdef NumberFieldElement x
        x = PY_NEW(NumberFieldElement)
        x._parent = <ParentWithBase>new_parent
        x.__denominator = self.__denominator
        cdef ntl_c_ZZX result
        cdef ntl_c_ZZ tmp
        cdef int i
        cdef int rel = _rel
        cdef ntl_ZZX _num
        cdef ntl_ZZ _den
        _num, _den = new_parent.polynomial_ntl()
        for i from 0 <= i <= deg(self.__numerator):
            GetCoeff(tmp, self.__numerator, i)
            SetCoeff(result, i*rel, tmp)
        rem_ZZX(x.__numerator, result, _num.x[0])
        return x

    def __reduce__(self):
        return __create__NumberFieldElement_version0, \
               (self.parent(), self.polynomial())

    def __repr__(self):
        x = self.polynomial()
        return str(x).replace(x.parent().variable_name(),self.parent().variable_name())

    def _im_gens_(self, codomain, im_gens):
        # NOTE -- if you ever want to change this so relative number fields are
        # in terms of a root of a poly.
        # The issue is that elements of a relative number field are represented in terms
        # of a generator for the absolute field.  However the morphism gives the image
        # of gen, which need not be a generator for the absolute field.  The morphism
        # has to be *over* the relative element.
        return codomain(self.polynomial()(im_gens[0]))

    def _latex_(self):
        """
        Returns the latex representation for this element.

        EXAMPLES:
            sage: C,zeta12=CyclotomicField(12).objgen()
            sage: latex(zeta12^4-zeta12)
            \zeta_{12}^{2} - \zeta_{12} - 1
        """
        return self.polynomial()._latex_(name=self.parent().latex_variable_name())

    def _pari_(self, var=None):
        """
        Return PARI C-library object corresponding to self.

        NOTE: This is not the actual underlying object that represents
        this element, since that is a polynomial in x (as PARI
        polynomials are rather constrained in their possible variable
        names, e.g., I cannot be the name of a variable).

        EXAMPLES:
            sage: k.<j> = QuadraticField(-1)
            sage: j._pari_()
            Mod(j, j^2 + 1)
            sage: pari(j)
            Mod(j, j^2 + 1)

        If you try do coerce a generator called I to PARI, hell may
        break loose:
            sage: k.<I> = QuadraticField(-1)
            sage: pari(I)
            Traceback (most recent call last):
            ...
            PariError: forbidden (45)

        Instead, request the variable be named different for the coercion:
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
            var = self.parent().variable_name()
        if isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_extension):
            f = self.polynomial()._pari_()
            g = str(self.parent().pari_relative_polynomial())
            base = self.parent().base_ring()
            gsub = base.gen()._pari_()
            gsub = str(gsub).replace(base.variable_name(), "y")
            g = g.replace("y", gsub)
        else:
            f = self.polynomial()._pari_()
            g = self.parent().polynomial()._pari_()
            if var != 'x':
                f = f.subst("x",var)
                g = g.subst("x",var)
        h = f.Mod(g)
        self.__pari[var] = h
        return h

    def _pari_init_(self, var=None):
        """
        Return GP/PARI string representation of self.
        """
        if var == None:
            var = self.parent().variable_name()
        if isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_extension):
            f = self.polynomial()._pari_()
            g = str(self.parent().pari_relative_polynomial())
            base = self.parent().base_ring()
            gsub = base.gen()._pari_()
            gsub = str(gsub).replace(base.variable_name(), "y")
            g = g.replace("y", gsub)
        else:
            f = self.polynomial()._pari_().subst("x",var)
            g = self.parent().polynomial()._pari_().subst("x",var)
        return 'Mod(%s, %s)'%(f,g)

    def __getitem__(self, n):
        return self.polynomial()[n]

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef NumberFieldElement _right = right
        return not (ZZX_equal(&left.__numerator, &_right.__numerator) and ZZ_equal(&left.__denominator, &_right.__denominator))

    def __abs__(self):
        return self.abs(i=0, prec=53)

    def abs(self, i=0, prec=53):
        """
        Return the absolute value of this element with respect to the
        ith complex embedding of parent, to the given precision.

        EXAMPLES:
            sage: z = CyclotomicField(7).gen()
            sage: abs(z)
            1.00000000000000
            sage: abs(z^2 + 17*z - 3)
            16.0604426799931
			sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3+17)
            sage: abs(a)
            2.57128159065824
            sage: a.abs(prec=100)
            2.5712815906582353554531872087
            sage: a.abs(1,100)
            2.5712815906582353554531872087
            sage: a.abs(2,100)
            2.5712815906582353554531872087

        Here's one where the absolute value depends on the embedding.
            sage: K.<b> = NumberField(x^2-2)
            sage: a = 1 + b
            sage: a.abs(i=0)
            2.41421356237309
            sage: a.abs(i=1)
            0.414213562373095
        """
        P = self.parent().complex_embeddings(prec)[i]
        return abs(P(self))

    def complex_embeddings(self, prec=53):
        phi = self.parent().complex_embeddings(prec)
        return [f(self) for f in phi]

    def complex_embedding(self, prec=53, i=0):
        return self.parent().complex_embeddings(prec)[i](self)

    def __pow__(self, r, mod):
        right = int(r)
        if right != r:
            raise NotImplementedError, "number field element to non-integral power not implemented"
        if right < 0:
            x = self.__invert__()
            right *= -1
            return sage.rings.arith.generic_power(x, right, one=self.parent()(1))
        return sage.rings.arith.generic_power(self, right, one=self.parent()(1))

    cdef void _reduce_c_(self):
        """
            Pull out common factors from the numerator and denominator!
        """
        cdef ntl_c_ZZ gcd
        cdef ntl_c_ZZ t1
        cdef ntl_c_ZZX t2
        content(t1, self.__numerator)
        GCD_ZZ(gcd, t1, self.__denominator)
        if sign(gcd) != sign(self.__denominator):
            negate(t1, gcd)
            gcd = t1
        div_ZZX_ZZ(t2, self.__numerator, gcd)
        div_ZZ_ZZ(t1, self.__denominator, gcd)
        self.__numerator = t2
        self.__denominator = t1

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        x = self._new()
        mul_ZZ(x.__denominator, self.__denominator, _right.__denominator)
        cdef ntl_c_ZZX t1, t2
        mul_ZZX_ZZ(t1, self.__numerator, _right.__denominator)
        mul_ZZX_ZZ(t2, _right.__numerator, self.__denominator)
        add_ZZX(x.__numerator, t1, t2)
        x._reduce_c_()
        return x

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        x = self._new()
        mul_ZZ(x.__denominator, self.__denominator, _right.__denominator)
        cdef ntl_c_ZZX t1, t2
        mul_ZZX_ZZ(t1, self.__numerator, _right.__denominator)
        mul_ZZX_ZZ(t2, _right.__numerator, self.__denominator)
        sub_ZZX(x.__numerator, t1, t2)
        x._reduce_c_()
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Returns the product of self and other as elements of a number field.

        EXAMPLES:
            sage: C.<zeta12>=CyclotomicField(12)
            sage: zeta12*zeta12^11
            1
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        x = self._new()
        mul_ZZ(x.__denominator, self.__denominator, _right.__denominator)
        cdef ntl_c_ZZ parent_den
        cdef ntl_c_ZZX parent_num
        self._parent_poly_c_( &parent_num, &parent_den )
        MulMod_ZZX(x.__numerator, self.__numerator, _right.__numerator, parent_num)
        x._reduce_c_()
        return x

        #NOTES: In LiDIA, they build a multiplication table for the
        #number field, so it's not necessary to reduce modulo the
        #defining polynomial every time:
        #     src/number_fields/algebraic_num/order.cc: compute_table
        # but asymptotically fast poly multiplication means it's
        # actually faster to *not* build a table!?!

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Returns the product of self and other as elements of a number field.

        EXAMPLES:
            sage: C.<I>=CyclotomicField(4)
            sage: 1/I
            -I
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        cdef ntl_c_ZZX inv_num
        cdef ntl_c_ZZ inv_den
        _right._invert_c_(&inv_num, &inv_den)
        x = self._new()
        mul_ZZ(x.__denominator, self.__denominator, inv_den)
        cdef ntl_c_ZZ parent_den
        cdef ntl_c_ZZX parent_num
        self._parent_poly_c_( &parent_num, &parent_den )
        MulMod_ZZX(x.__numerator, self.__numerator, inv_num, parent_num)
        x._reduce_c_()
        return x

    def __floordiv__(self, other):
        return self / other

    cdef ModuleElement _neg_c_impl(self):
        cdef NumberFieldElement x
        x = self._new()
        mul_ZZX_long(x.__numerator, self.__numerator, -1)
        x.__denominator = self.__denominator
        return x

    def __int__(self):
        """
        EXAMPLES:
            sage: C.<I>=CyclotomicField(4)
            sage: int(1/I)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
            sage: int(I*I)
            -1
        """
        return int(self.polynomial())

    def __long__(self):
        return long(self.polynomial())

    cdef void _parent_poly_c_(self, ntl_c_ZZX *num, ntl_c_ZZ *den):
        cdef long i
        cdef ntl_c_ZZ coeff
        cdef ntl_ZZX _num
        cdef ntl_ZZ _den
        if isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_extension):
            # ugly temp code
            f = self.parent().absolute_polynomial()

            __den = f.denominator()
            (<Integer>ZZ(__den))._to_ZZ(den)

            __num = f * __den
            for i from 0 <= i <= __num.degree():
                (<Integer>ZZ(__num[i]))._to_ZZ(&coeff)
                SetCoeff( num[0], i, coeff )
        else:
            _num, _den = self.parent().polynomial_ntl()
            num[0] = _num.x[0]
            den[0] = _den.x[0]

    cdef void _invert_c_(self, ntl_c_ZZX *num, ntl_c_ZZ *den):
        """
        Computes the numerator and denominator of the multiplicative inverse of this element.

        Suppose that this element is x/d and the parent mod'ding polynomial is M/D.  The NTL function
        XGCD( r, s, t, a, b ) computes r,s,t such that $r=s*a+t*b$.  We compute
        XGCD( r, s, t, x*D, M*d ) and set
        num=s*D*d
        den=r

        EXAMPLES:
            I'd love to, but since we are dealing with c-types, I can't at this level.
            Check __invert__ for doc-tests that rely on this functionality.
        """
        cdef ntl_c_ZZ parent_den
        cdef ntl_c_ZZX parent_num
        self._parent_poly_c_( &parent_num, &parent_den )

        cdef ntl_c_ZZX t # unneeded except to be there
        cdef ntl_c_ZZX a, b
        mul_ZZX_ZZ( a, self.__numerator, parent_den )
        mul_ZZX_ZZ( b, parent_num, self.__denominator )
        XGCD_ZZX( den[0], num[0],  t, a, b, 1 )
        mul_ZZX_ZZ( num[0], num[0], parent_den )
        mul_ZZX_ZZ( num[0], num[0], self.__denominator )

    def __invert__(self):
        """
        Returns the multiplicative inverse of self in the number field.

        EXAMPLES:
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
#        if isinstance(K, sage.rings.number_field.number_field.NumberField_extension):
#            return K(K.pari_rnf().rnfeltreltoabs(quotient))
#        else:
#            return K(quotient)

    def _integer_(self):
        """
        Returns a rational integer if this element is actually a rational integer.

        EXAMPLES:
            sage: C.<I>=CyclotomicField(4)
            sage: (~I)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce -I to an integer
            sage: (2*I*I)._integer_()
            -2
        """
        if deg(self.__numerator) >= 1:
            raise TypeError, "Unable to coerce %s to an integer"%self
        return ZZ(self._rational_())

    def _rational_(self):
        """
        Returns a rational number if this element is actually a rational number.

        EXAMPLES:
            sage: C.<I>=CyclotomicField(4)
            sage: (~I)._rational_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce -I to a rational
            sage: (I*I/2)._rational_()
            -1/2
        """
        if deg(self.__numerator) >= 1:
            raise TypeError, "Unable to coerce %s to a rational"%self
        cdef Integer num
        num = PY_NEW(Integer)
        ZZX_getitem_as_mpz(&num.value, &self.__numerator, 0)
        return num / (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)

    def conjugate(self):
        """
        Return the complex conjugate of the number field element.  Currently,
        this is implemented for cyclotomic fields and quadratic extensions of Q.
        It seems likely that there are other number fields for which the idea of
        a conjugate would be easy to compute.

        EXAMPLES:
            sage: k.<I> = QuadraticField(-1)
            sage: I.conjugate()
            -I
            sage: (I/(1+I)).conjugate()
            -1/2*I + 1/2
            sage: z6=CyclotomicField(6).gen(0)
            sage: (2*z6).conjugate()
            -2*zeta6 + 2

			sage: x = polygen(QQ, 'x')
            sage: K.<b> = NumberField(x^3 - 2)
            sage: b.conjugate()
            Traceback (most recent call last):
            ...
            NotImplementedError: complex conjugation is not implemented (or doesn't make sense).
        """
        coeffs = self.parent().polynomial().list()
        if len(coeffs) == 3 and coeffs[2] == 1 and coeffs[1] == 0:
            # polynomial looks like x^2+d
            # i.e. we live in a quadratic extension of QQ
            if coeffs[0] > 0:
                gen = self.parent().gen()
                return self.polynomial()(-gen)
            else:
                return self
        elif isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_cyclotomic):
            # We are in a cyclotomic field
            # Replace the generator zeta_n with (zeta_n)^(n-1)
            gen = self.parent().gen()
            return self.polynomial()(gen ** (gen.multiplicative_order()-1))
        else:
            raise NotImplementedError, "complex conjugation is not implemented (or doesn't make sense)."

    def polynomial(self):
        coeffs = []
        cdef Integer den = (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)
        cdef Integer numCoeff
        cdef int i
        for i from 0 <= i <= deg(self.__numerator):
            numCoeff = PY_NEW(Integer)
            ZZX_getitem_as_mpz(&numCoeff.value, &self.__numerator, i)
            coeffs.append( numCoeff / den )
        return QQ['x'](coeffs)

    def denominator(self):
        """
        Return the denominator of this element, which is by definition
        the denominator of the corresponding polynomial
        representation.  I.e., elements of number fields are
        represented as a polynomial (in reduced form) modulo the
        modulus of the number field, and the denominator is the
        denominator of this polynomial.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: a = 1/3 + (1/5)*z
            sage: print a.denominator()
            15
        """
        return (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)

    def _set_multiplicative_order(self, n):
        self.__multiplicative_order = n

    def multiplicative_order(self):
        if self.__multiplicative_order is not None:
            return self.__multiplicative_order

        if deg(self.__numerator) == 0:
            if self._rational_() == 1:
                self.__multiplicative_order = 1
                return self.__multiplicative_order
            if self._rational_() == -1:
                self.__multiplicative_order = 2
                return self.__multiplicative_order

        if isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_cyclotomic):
            t = self.parent().multiplicative_order_table()
            f = self.polynomial()
            if t.has_key(f):
                self.__multiplicative_order = t[f]
                return self.__multiplicative_order

        ####################################################################
        # VERY DUMB Algorithm to compute the multiplicative_order of
        # an element x of a number field K.
        #
        # 1. Find an integer B such that if n>=B then phi(n) > deg(K).
        #    For this use that for n>6 we have phi(n) >= log_2(n)
        #    (to see this think about the worst prime factorization
        #    in the multiplicative formula for phi.)
        # 2. Compute x, x^2, ..., x^B in order to determine the multiplicative_order.
        #
        # todo-- Alternative: Only do the above if we don't require an optional
        # argument which gives a multiple of the order, which is usually
        # something available in any actual application.
        #
        # BETTER TODO: Factor cyclotomic polynomials over K to determine
        # possible orders of elements?  Is there something even better?
        #
        ####################################################################
        d = self.parent().degree()
        B = max(7, 2**d+1)
        x = self
        i = 1
        while i < B:
            if x == 1:
                self.__multiplicative_order = i
                return self.__multiplicative_order
            x *= self
            i += 1

        # it must have infinite order
        self.__multiplicative_order = infinity.infinity
        return self.__multiplicative_order

    def trace(self):
        K = self.parent().base_ring()
        return K(self._pari_('x').trace())

    def norm(self):
        K = self.parent().base_ring()
        return K(self._pari_('x').norm())

    def charpoly(self, var):
        r"""
        The characteristic polynomial of this element over $\Q$.

        EXAMPLES:

        We compute the charpoly of cube root of $3$.

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: a.charpoly('x')
            x^3 - 2

        We construct a relative extension and find the characteristic
        polynomial over $\Q$.

            sage: S.<X> = K[]
            sage: L.<b> = NumberField(X^3 + 17)
            sage: L
            Extension by X^3 + 17 of the Number Field in a with defining polynomial x^3 - 2
            sage: a = L.0; a
            b
            sage: a.charpoly('x')
            x^9 + 57*x^6 + 165*x^3 + 6859
            sage: a.charpoly('y')
            y^9 + 57*y^6 + 165*y^3 + 6859
        """
        R = self.parent().base_ring()[var]
        if not isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_extension):
            return R(self._pari_('x').charpoly())
        else:
            g = self.polynomial()  # in QQ[x]
            f = self.parent().pari_polynomial()  # # field is QQ[x]/(f)
            return R( (g._pari_('x').Mod(f)).charpoly() )

## This might be useful for computing relative charpoly.
## BUT -- currently I don't even know how to view elements
## as being in terms of the right thing, i.e., this code
## below as is lies.
##             nf = self.parent()._pari_base_nf()
##             prp = self.parent().pari_relative_polynomial()
##             elt = str(self.polynomial()._pari_())
##             return R(nf.rnfcharpoly(prp, elt))
##         # return self.matrix().charpoly('x')

    def minpoly(self, var='x'):
        """
        Return the minimal polynomial of this number field element.

        EXAMPLES:
			sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^2+3)
            sage: a.minpoly('x')
            x^2 + 3
            sage: R.<X> = K['X']
            sage: L.<b> = K.extension(X^2-(22 + a))
            sage: b.minpoly('t')
            t^4 + (-44)*t^2 + 487
            sage: b^2 - (22+a)
            0
        """
        # The minimal polynomial is square-free and
        # divisible by same irreducible factors as
        # the characteristic polynomial.
        # TODO: factoring to find the square-free part is idiotic.
        # Instead use a GCD algorithm!
        f = sage.rings.polynomial_ring.PolynomialRing(QQ, str(var))(1)
        for g, _ in self.charpoly(var).factor():
            f *= g
        return f

    def matrix(self):
        r"""
        The matrix of right multiplication by the element on the power
        basis $1, x, x^2, \ldots, x^{d-1}$ for the number field.  Thus
        the {\em rows} of this matrix give the images of each of the $x^i$.

        EXAMPLES:

        Regular number field:
            sage: K.<a> = NumberField(QQ['x'].0^3 - 5)
            sage: M = a.matrix(); M
            [0 1 0]
            [0 0 1]
            [5 0 0]
            sage: M.base_ring() is QQ
            True

        """
##         Relative number field:
##             sage: L.<b> = K.extension(K['x'].0^2 - 2)
##             sage: 1*b, b*b, b**3, b**6
##             (b, b^2, b^3, 6*b^4 - 10*b^3 - 12*b^2 - 60*b - 17)
##             sage: L.pari_rnf().rnfeltabstorel(b._pari_())
##             x - y
##             sage: L.pari_rnf().rnfeltabstorel((b**2)._pari_())
##             2
##             sage: M = b.matrix(); M
##             [0 1]
##             [3 0]
##             sage: M.base_ring() is K
##             True

#         Absolute number field:
#             sage: M = L.absolute_field().gen().matrix(); M
#             [  0   1   0   0   0   0]
#             [  0   0   1   0   0   0]
#             [  0   0   0   1   0   0]
#             [  0   0   0   0   1   0]
#             [  0   0   0   0   0   1]
#             [  2 -90 -27 -10   9   0]
#             sage: M.base_ring() is QQ
#             True

#         More complicated relative number field:
#             sage: L.<b> = K.extension(K['x'].0^2 - a); L
#             Extension by x^2 + -a of the Number Field in a with defining polynomial x^3 - 5
#             sage: M = b.matrix(); M
#             [0 1]
#             [a 0]
#             sage: M.base_ring()
#             sage: M.base_ring() is K
#             True
        # Mutiply each power of field generator on
        # the left by this element; make matrix
        # whose rows are the coefficients of the result,
        # and transpose.
        if self.__matrix is None:
            K = self.parent()
            v = []
            x = K.gen()
            a = K(1)
            d = K.degree()
            for n in range(d):
                v += (a*self).list()
                a *= x
            k = K.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(k, d)
            self.__matrix = M(v)
        return self.__matrix

    def list(self):
        """
        EXAMPLE:
            sage: K.<z> = CyclotomicField(3)
            sage: (2+3/5*z).list()
            [2, 3/5]
            sage: (5*z).list()
            [0, 5]
            sage: K(3).list()
            [3, 0]
        """
        P = self.parent()
        # The algorithm below is total nonsense, unless the parent of self is an
        # absolute extension.
        if isinstance(P, sage.rings.number_field.number_field.NumberField_extension):
            raise NotImplementedError
        n = self.parent().degree()
        v = self.polynomial().list()[:n]
        z = sage.rings.rational.Rational(0)
        return v + [z]*(n - len(v))
