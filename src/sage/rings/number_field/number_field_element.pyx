"""
Number Field Elements

AUTHORS:
    -- William Stein version before it got cython'd
    -- Joel B. Mohler (2007-03-09): First reimplementation into cython
    -- William Stein (2007-09-04): add doctests
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

include "../../ext/stdsage.pxi"
include '../../ext/interrupt.pxi'

import sage.rings.field_element
import sage.rings.infinity
import sage.rings.polynomial.polynomial_element
import sage.rings.polynomial.polynomial_ring
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
    """
    Return True if x is of type NumberFieldElement, i.e., an
    element of a number field.

    EXAMPLES:
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

        sage: k.<a> = NumberField(x^3 - 2)
        sage: loads(dumps(a+1)) == a + 1
        True
    """
    return NumberFieldElement(parent, poly)

cdef class NumberFieldElement(FieldElement):
    """
    An element of a number field.

    EXAMPLES:
        sage: k.<a> = NumberField(x^3 + x + 1)
        sage: a^3
        -a - 1
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

        cdef ZZ_c coeff
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
        if not isinstance(f, sage.rings.polynomial.polynomial_element.Polynomial):
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

    def __alloc__(self):
        ZZX_construct(&self.__numerator)
        ZZ_construct(&self.__denominator)

    def __dealloc__(self):
        ZZX_destruct(&self.__numerator)
        ZZ_destruct(&self.__denominator)

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
        cdef ZZX_c result
        cdef ZZ_c tmp
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
        """
        Used in pickling number field elements.

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 - 17*x^2 + 1)
            sage: t = a.__reduce__(); t
            (<built-in function __create__NumberFieldElement_version0>, (Number Field in a with defining polynomial x^3 - 17*x^2 + 1, x))
            sage: t[0](*t[1]) == a
            True
        """
        return __create__NumberFieldElement_version0, \
               (self.parent(), self.polynomial())

    def __repr__(self):
        """
        String representation of this number field element,
        which is just a polynomial in the generator.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 2)
            sage: b = (2/3)*a + 3/5
            sage: b.__repr__()
            '2/3*a + 3/5'
        """
        x = self.polynomial()
        return str(x).replace(x.parent().variable_name(),self.parent().variable_name())

    def _im_gens_(self, codomain, im_gens):
        """
        This is used in computing homomorphisms between number fields.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 - 2)
            sage: m.<b> = NumberField(x^4 - 2)
            sage: phi = k.hom([b^2])
            sage: phi(a+1)
            b^2 + 1
            sage: (a+1)._im_gens_(m, [b^2])
            b^2 + 1
        """
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

    def _pari_(self, var='x'):
        """
        Return PARI C-library object corresponding to self.

        EXAMPLES:
            sage: k.<j> = QuadraticField(-1)
            sage: j._pari_('j')
            Mod(j, j^2 + 1)
            sage: pari(j)
            Mod(x, x^2 + 1)

            sage: y = QQ['y'].gen()
            sage: k.<j> = NumberField(y^3 - 2)
            sage: pari(j)
            Mod(x, x^3 - 2)

        By default the variable name is 'x', since in PARI many variable
        names are reserved:
            sage: theta = polygen(QQ, 'theta')
            sage: M.<theta> = NumberField(theta^2 + 1)
            sage: pari(theta)
            Mod(x, x^2 + 1)

        If you try do coerce a generator called I to PARI, hell may
        break loose:
            sage: k.<I> = QuadraticField(-1)
            sage: I._pari_('I')
            Traceback (most recent call last):
            ...
            PariError: forbidden (45)

        Instead, request the variable be named different for the coercion:
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
            var = self.parent().variable_name()
        if isinstance(self.parent(),
                      sage.rings.number_field.number_field.NumberField_extension):
            f = self.polynomial()._pari_()
            g = str(self.parent().pari_polynomial())
            base = self.parent().base_ring()
            gsub = base.gen()._pari_()
            gsub = str(gsub).replace('x', "y")
            g = g.replace("y", gsub)
        else:
            f = self.polynomial()._pari_()
            gp = self.parent().polynomial()
            if gp.name() != 'x':
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

    def _pari_init_(self, var='x'):
        """
        Return GP/PARI string representation of self. This is used for
        converting this number field element to GP/PARI.  The returned
        string defines a pari Mod in the variable is var, which is by
        default 'x' -- not the name of the generator of the number
        field.

        INPUT:
            var -- (default: 'x') the variable of the pari Mod.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 - x - 1)
            sage: ((1 + 1/3*a)^4)._pari_init_()
            'Mod(1/81*x^4 + 4/27*x^3 + 2/3*x^2 + 4/3*x + 1, x^5 - x - 1)'
            sage: ((1 + 1/3*a)^4)._pari_init_('a')
            'Mod(1/81*a^4 + 4/27*a^3 + 2/3*a^2 + 4/3*a + 1, a^5 - a - 1)'

        Note that _pari_init_ can fail because of reserved words in PARI,
        and since it actually works by obtaining the PARI representation
        of something.
            sage: K.<theta> = NumberField(x^5 - x - 1)
            sage: b = (1/2 - 2/3*theta)^3; b
            -8/27*theta^3 + 2/3*theta^2 - 1/2*theta + 1/8
            sage: b._pari_init_('theta')
            Traceback (most recent call last):
            ...
            PariError: unexpected character (2)

        Fortunately pari_init returns everything in terms of x by default.
            sage: pari(b)
            Mod(-8/27*x^3 + 2/3*x^2 - 1/2*x + 1/8, x^5 - x - 1)
        """
        return repr(self._pari_(var=var))
##         if var == None:
##             var = self.parent().variable_name()
##         if isinstance(self.parent(), sage.rings.number_field.number_field.NumberField_extension):
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

        Note that $n$ must be between 0 and $d-1$, where $d$ is the
        degree of the number field.

        EXAMPLES:
            sage: m.<b> = NumberField(x^4 - 1789)
            sage: c = (2/3-4/5*b)^3; c
            -64/125*b^3 + 32/25*b^2 - 16/15*b + 8/27
            sage: c[0]
            8/27
            sage: c[2]
            32/25
            sage: c[3]
            -64/125

        We illustrate bounds checking:
            sage: c[-1]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1.
            sage: c[4]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1.

        The list method implicitly calls __getitem__:
            sage: list(c)
            [8/27, -16/15, 32/25, -64/125]
            sage: m(list(c)) == c
            True
        """
        if n < 0 or n >= self.parent().degree():     # make this faster.
            raise IndexError, "index must be between 0 and degree minus 1."
        return self.polynomial()[n]

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef NumberFieldElement _right = right
        return not (ZZX_equal(&left.__numerator, &_right.__numerator) and ZZ_equal(&left.__denominator, &_right.__denominator))

    def __abs__(self):
        r"""
        Return the numerical absolute value of this number field
        element with respect to the first archimedean embedding, to 53
        bits of precision.

        This is the \code{abs( )} Python function.  If you want a different
        embedding or precision, use \code{self.abs(...)}.

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 - 2)
            sage: abs(a)
            1.25992104989487
            sage: abs(a)^3
            2.00000000000000
            sage: a.abs(prec=128)
            1.2599210498948731647672106072782283506
        """
        return self.abs(prec=53, i=0)

    def abs(self, prec=53, i=0):
        """
        Return the absolute value of this element with respect to the
        ith complex embedding of parent, to the given precision.

        INPUT:
            prec -- (default: 53) integer bits of precision
            i -- (default: ) integer, which embedding to use

        EXAMPLES:
            sage: z = CyclotomicField(7).gen()
            sage: abs(z)
            1.00000000000000
            sage: abs(z^2 + 17*z - 3)
            16.0604426799931
            sage: K.<a> = NumberField(x^3+17)
            sage: abs(a)
            2.57128159065824
            sage: a.abs(prec=100)
            2.5712815906582353554531872087
            sage: a.abs(prec=100,i=1)
            2.5712815906582353554531872087
            sage: a.abs(100, 2)
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
        """
        Return the images of this element in the floating point
        complex numbers, to the given bits of precision.

        INPUT:
            prec -- integer (default: 53) bits of precision

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 - 2)
            sage: a.complex_embeddings()
            [1.25992104989487, -0.629960524947437 + 1.09112363597172*I, -0.629960524947437 - 1.09112363597172*I]
            sage: a.complex_embeddings(10)
            [1.3, -0.63 + 1.1*I, -0.63 - 1.1*I]
            sage: a.complex_embeddings(100)
            [1.2599210498948731647672106073, -0.62996052494743658238360530364 + 1.0911236359717214035600726142*I, -0.62996052494743658238360530364 - 1.0911236359717214035600726142*I]
        """
        phi = self.parent().complex_embeddings(prec)
        return [f(self) for f in phi]

    def complex_embedding(self, prec=53, i=0):
        """
        Return the i-th embedding of self in the complex numbers, to
        the given precision.

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 - 2)
            sage: a.complex_embedding()
            1.25992104989487
            sage: a.complex_embedding(10)
            1.3
            sage: a.complex_embedding(100)
            1.2599210498948731647672106073
            sage: a.complex_embedding(20, 1)
            -0.62996 + 1.0911*I
            sage: a.complex_embedding(20, 2)
            -0.62996 - 1.0911*I
        """
        return self.parent().complex_embeddings(prec)[i](self)

    def is_square(self, root=False):
        """
        Return True if self is a square in its parent number field and
        otherwise return False.

        INPUT:
            root -- if True, also return a square root (or None if self
                    is not a perfect square)

        EXAMPLES:
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

        EXAMPLES:
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

            sage: K.<a> = CyclotomicField(7)
            sage: a.sqrt()
            a^4

            sage: K.<a> = NumberField(x^5 - x + 1)
            sage: (a^4 + a^2 - 3*a + 2).sqrt()
            a^3 - a^2

        ALGORITHM:
            Use Pari to factor $x^2$ - \code{self} in K.

        """
        # For now, use pari's factoring abilities
        R = sage.rings.polynomial.polynomial_ring.PolynomialRing(self._parent, 't')
        f = R([-self, 0, 1])
        roots = f.roots()
        if all:
            return [r[0] for r in roots]
        elif len(roots) > 0:
            return roots[0][0]
        else:
            raise ValueError, "%s not a square in %s"%(self, self._parent)

    cdef void _reduce_c_(self):
        """
        Pull out common factors from the numerator and denominator!
        """
        cdef ZZ_c gcd
        cdef ZZ_c t1
        cdef ZZX_c t2
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
        cdef ZZX_c t1, t2
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
        cdef ZZX_c t1, t2
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
        cdef ZZ_c parent_den
        cdef ZZX_c parent_num
        self._parent_poly_c_( &parent_num, &parent_den )
        x = self._new()
        _sig_on
        # MulMod doesn't handle non-monic polynomials.
        # Therefore, we handle the non-monic case entirely separately.
        if ZZX_is_monic( &parent_num ):
            mul_ZZ(x.__denominator, self.__denominator, _right.__denominator)
            MulMod_ZZX(x.__numerator, self.__numerator, _right.__numerator, parent_num)
        else:
            mul_ZZ(x.__denominator, self.__denominator, _right.__denominator)
            mul_ZZX(x.__numerator, self.__numerator, _right.__numerator)
            if ZZX_degree(&x.__numerator) >= ZZX_degree(&parent_num):
                mul_ZZX_ZZ( x.__numerator, x.__numerator, parent_den )
                mul_ZZX_ZZ( temp, parent_num, x.__denominator )
                power_ZZ(temp1,LeadCoeff_ZZX(temp),ZZX_degree(&x.__numerator)-ZZX_degree(&parent_num)+1)
                PseudoRem_ZZX(x.__numerator, x.__numerator, temp)
                mul_ZZ(x.__denominator, x.__denominator, parent_den)
                mul_ZZ(x.__denominator, x.__denominator, temp1)
        _sig_off
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
        Returns the quotient of self and other as elements of a number field.

        EXAMPLES:
            sage: C.<I>=CyclotomicField(4)
            sage: 1/I
            -I
            sage: I/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Number field element division by zero

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
        cdef ZZ_c parent_den
        cdef ZZX_c parent_num
        cdef ZZX_c temp
        cdef ZZ_c temp1
        if not _right:
            raise ZeroDivisionError, "Number field element division by zero"
        self._parent_poly_c_( &parent_num, &parent_den )
        x = self._new()
        _sig_on
        _right._invert_c_(&inv_num, &inv_den)
        if ZZX_is_monic( &parent_num ):
            mul_ZZ(x.__denominator, self.__denominator, inv_den)
            MulMod_ZZX(x.__numerator, self.__numerator, inv_num, parent_num)
        else:
            mul_ZZ(x.__denominator, self.__denominator, inv_den)
            mul_ZZX(x.__numerator, self.__numerator, inv_num)
            if ZZX_degree(&x.__numerator) >= ZZX_degree(&parent_num):
                mul_ZZX_ZZ( x.__numerator, x.__numerator, parent_den )
                mul_ZZX_ZZ( temp, parent_num, x.__denominator )
                power_ZZ(temp1,LeadCoeff_ZZX(temp),ZZX_degree(&x.__numerator)-ZZX_degree(&parent_num)+1)
                PseudoRem_ZZX(x.__numerator, x.__numerator, temp)
                mul_ZZ(x.__denominator, x.__denominator, parent_den)
                mul_ZZ(x.__denominator, x.__denominator, temp1)
        x._reduce_c_()
        _sig_off
        return x

    def __floordiv__(self, other):
        """
        Return the quotient of self and other.  Since these are field
        elements the floor division is exactly the same as usual
        division.

        EXAMPLES:
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

        EXAMPLES:
            sage: m.<b> = CyclotomicField(17)
            sage: m(0).__nonzero__()
            False
            sage: b.__nonzero__()
            True

        Nonzero is used by the bool command:
            sage: bool(b + 1)
            True
            sage: bool(m(0))
            False
        """
        return not IsZero_ZZX(self.__numerator)

    cdef ModuleElement _neg_c_impl(self):
        cdef NumberFieldElement x
        x = self._new()
        mul_ZZX_long(x.__numerator, self.__numerator, -1)
        x.__denominator = self.__denominator
        return x

    def __int__(self):
        """
        Attempt to convert this number field element to a Python integer,
        if possible.

        EXAMPLES:
            sage: C.<I>=CyclotomicField(4)
            sage: int(1/I)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
            sage: int(I*I)
            -1

            sage: K.<a> = NumberField(x^10 - x - 1)
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
            sage: int(K(9390283))
            9390283

        The semantics are like in Python, so the value does not have
        to preserved.
            sage: int(K(393/29))
            13
        """
        return int(self.polynomial())

    def __long__(self):
        """
        Attempt to convert this number field element to a Python long,
        if possible.

        EXAMPLES:
            sage: K.<a> = NumberField(x^10 - x - 1)
            sage: long(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to long
            sage: long(K(1234))
            1234L

        The value does not have to be preserved, in the case of fractions.
            sage: long(K(393/29))
            13L
        """
        return long(self.polynomial())

    cdef void _parent_poly_c_(self, ZZX_c *num, ZZ_c *den):
        cdef long i
        cdef ZZ_c coeff
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

    cdef void _invert_c_(self, ZZX_c *num, ZZ_c *den):
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
        cdef ZZ_c parent_den
        cdef ZZX_c parent_num
        self._parent_poly_c_( &parent_num, &parent_den )

        cdef ZZX_c t # unneeded except to be there
        cdef ZZX_c a, b
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
        """
        Return the underlyling polynomial corresponding to this
        number field element.

        The resulting polynomial is currently *not* cached.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 - x - 1)
            sage: f = (-2/3 + 1/3*a)^4; f
            1/81*a^4 - 8/81*a^3 + 8/27*a^2 - 32/81*a + 16/81
            sage: g = f.polynomial(); g
            1/81*x^4 - 8/81*x^3 + 8/27*x^2 - 32/81*x + 16/81
            sage: parent(g)
            Univariate Polynomial Ring in x over Rational Field

        Note that the result of this function is not cached (should this
        be changed?):
            sage: g is f.polynomial()
            False
        """
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
        """
        Set the multiplicative order of this number field element.

        WARNING -- use with caution -- only for internal use!  End
        users should never call this unless they have a very good
        reason to do so.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: a._set_multiplicative_order(3)
            sage: a.multiplicative_order()
            3

        You can be evil with this so be careful.  That's why the function
        name begins with an underscore.
            sage: a._set_multiplicative_order(389)
            sage: a.multiplicative_order()
            389
        """
        self.__multiplicative_order = n

    def multiplicative_order(self):
        """
        Return the multiplicative order of this number field element.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(5)
            sage: z.multiplicative_order()
            5
            sage: (-z).multiplicative_order()
            10
            sage: (1+z).multiplicative_order()
            +Infinity
        """
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
            t = self.parent()._multiplicative_order_table()
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
        self.__multiplicative_order = sage.rings.infinity.infinity
        return self.__multiplicative_order

    def trace(self):
        """
        Return the trace of this number field element.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 -132/7*x^2 + x + 1); K
            Number Field in a with defining polynomial x^3 - 132/7*x^2 + x + 1
            sage: a.trace()
            132/7
            sage: (a+1).trace() == a.trace() + 3
            True
        """
        K = self.parent().base_ring()
        return K(self._pari_('x').trace())

    def norm(self):
        """
        Return the norm of this number field element.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + x^2 + x + -132/7); K
            Number Field in a with defining polynomial x^3 + x^2 + x - 132/7
            sage: a.norm()
            132/7
            sage: K(0).norm()
            0
        """
        K = self.parent().base_ring()
        return K(self._pari_('x').norm())

    def charpoly(self, var='x'):
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
            sage: L.<b> = NumberField(X^3 + 17); L
            Number Field in b with defining polynomial X^3 + 17 over its base field
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
            return R( (g._pari_().Mod(f)).charpoly() )

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
        return self.charpoly(var).radical() # square free part of charpoly

    def is_integral(self):
        r"""
        Determine if a number is in the ring of integers
        of this number field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 23, 'a')
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
        """
        return all([a in ZZ for a in self.minpoly()])

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
