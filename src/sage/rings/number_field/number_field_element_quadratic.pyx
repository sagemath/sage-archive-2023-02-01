"""nodoctest
Not done.
"""
#*****************************************************************************
#     Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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

"""
Optimized Quadratic Number Field Elements

NO TESTS => BROKEN!

AUTHORS:
    -- Robert Bradshaw (2007-09): Initial version
"""

include '../../ext/interrupt.pxi'
include "../../ext/stdsage.pxi"
cdef extern from *:
    # TODO: move to stdsage.pxi
    object PY_NEW_SAME_TYPE(object o)
    bint PY_TYPE_CHECK_EXACT(object o, object type)

cdef extern from *:
    # TODO: move to cdefs.pxi
    void mpz_addmul (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_submul (mpz_t rop, mpz_t op1, mpz_t op2)

cdef object QQ, ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ

cdef mpz_to_str(mpz_t z):
    cdef Integer zz = PY_NEW(Integer)
    mpz_set(zz.value, z)
    return str(zz)

def __make_NumberFieldElement_quadratic0(parent, a, b, denom):
    """
    Used in unpickling elements of number fields.

    TEST:
        sage: K.<a> = NumberField(x^2-x+13)
        sage: loads(dumps(a)) == a
        True
    """
    return NumberFieldElement_quadratic(parent, (a, b, denom))


cdef class NumberFieldElement_quadratic(NumberFieldElement_absolute):
    # (a + b sqrt(disc)) / denom
    def __init__(self, parent, f):
        """
        Construct a NumberFieldElement_quadratic object as an efficiently
        represented member of an absolute quadratic field.

        TESTS:
            sage: from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic

          Setup some fields:
            sage: K.<a> = NumberField(x^2+23)
            sage: a.parts()
            (0, 1)
            sage: F.<b> = NumberField(x^2-x+7)
            sage: b.parts()
            (1/2, 3/2)

          By polynomials:
            sage: NumberFieldElement_quadratic(K, x-1)
            a - 1
            sage: NumberFieldElement_quadratic(F, x-1)
            b - 1

          By triples of Integers:
            sage: NumberFieldElement_quadratic(K, (1,2,3))
            2/3*a + 1/3
            sage: NumberFieldElement_quadratic(F, (1,2,3))
            4/9*b + 1/9
            sage: NumberFieldElement_quadratic(F, (1,2,3)).parts()
            (1/3, 2/3)

          By pairs of Rationals:
            sage: NumberFieldElement_quadratic(K, (1/2,1/3))
            1/3*a + 1/2
            sage: NumberFieldElement_quadratic(F, (1/2,1/3))
            2/9*b + 7/18
            sage: NumberFieldElement_quadratic(F, (1/2,1/3)).parts()
            (1/2, 1/3)

          Direct from Rational:
            sage: NumberFieldElement_quadratic(K, 2/3)
            2/3
            sage: NumberFieldElement_quadratic(F, 2/3)
            2/3
            """
        self.disc = parent.discriminant()
        cdef Integer a, b, denom
        cdef Rational ad, bd

        cdef NumberFieldElement_quadratic gen

        if PY_TYPE_CHECK(f, NumberFieldElement_quadratic):
            mpz_set(self.a, (<NumberFieldElement_quadratic>f).a)
            mpz_set(self.b, (<NumberFieldElement_quadratic>f).b)
            mpz_set(self.denom, (<NumberFieldElement_quadratic>f).denom)

        elif PY_TYPE_CHECK_EXACT(f, Rational):
            NumberFieldElement_absolute.__init__(self, parent, None)
            mpz_set(self.a, mpq_numref((<Rational>f).value))
            mpz_set_ui(self.b, 0)
            mpz_set(self.denom, mpq_denref((<Rational>f).value))

        elif PY_TYPE_CHECK_EXACT(f, tuple) and len(f) == 2:
            NumberFieldElement_absolute.__init__(self, parent, None)
            ad, bd = f
            mpz_set(self.a, a.value)
            mpz_set(self.b, b.value)
            mpz_lcm(self.denom, mpq_denref(ad.value), mpq_denref(bd.value))
            mpz_divexact(self.a, self.denom, mpq_denref(ad.value))
            mpz_mul(self.a, self.a, mpq_numref(ad.value))
            mpz_divexact(self.b, self.denom, mpq_denref(bd.value))
            mpz_mul(self.b, self.b, mpq_numref(bd.value))

        elif PY_TYPE_CHECK_EXACT(f, tuple) and len(f) == 3:
            NumberFieldElement_absolute.__init__(self, parent, None)
            a, b, denom = f
            mpz_set(self.a, a.value)
            mpz_set(self.b, b.value)
            mpz_set(self.denom, denom.value)
            self._reduce_c_()

        else:
            NumberFieldElement_absolute.__init__(self, parent, f)
            # poly is in gen (which may not be sqrt(d))
            self._ntl_coeff_as_mpz(&self.a, 0)
            self._ntl_coeff_as_mpz(&self.b, 1)
            if mpz_cmp_ui(self.a, 0) or mpz_cmp_ui(self.b, 0):
                gen = parent.gen() # should this be cached?
                self._ntl_denom_as_mpz(&self.denom)
                if mpz_cmp_ui(self.b, 0):
                    mpz_mul(self.a, self.a, gen.denom)
                    mpz_addmul(self.a, self.b, gen.a)
                    mpz_mul(self.b, self.b, gen.b)
                    mpz_mul(self.denom, self.denom, gen.denom)
            else:
                mpz_set_ui(self.denom, 1)

    def __new__(self, parent=None, f=None):
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.denom)

    def __dealloc__(self):
        mpz_clear(self.a)
        mpz_clear(self.b)
        mpz_clear(self.denom)

    def __reduce__(self):
        """
        TEST:
            sage: K.<a> = NumberField(x^2-13)
            sage: loads(dumps(a)) == a
            True
            sage: loads(dumps(a/3+5)) == a/3+5
            True
        """
        cdef Integer a = <Integer>PY_NEW(Integer)
        cdef Integer b = <Integer>PY_NEW(Integer)
        cdef Integer denom = <Integer>PY_NEW(Integer)
        mpz_set(a.value, self.a)
        mpz_set(b.value, self.b)
        mpz_set(denom.value, self.denom)
        return __make_NumberFieldElement_quadratic0, (self._parent, a, b, denom)

    def parts(self):
        """
        This function returns the internal representation of self as
        a pair of rationals $a$ and $b$, where self = $a+b\sqrt{D}$.

        This is much closer to the internal storage format of the
        elements, unless the generator is equal to $\sqrt{D}$ will
        be different than the polynomial representation coefficients.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2-13)
            sage: K.discriminant()
            13
            sage: a.parts()
            (0, 1)
            sage: (a/2-4).parts()
            (-4, 1/2)
            sage: K.<a> = NumberField(x^2-7)
            sage: K.discriminant()
            28
            sage: a.parts()
            (0, 1/2)
            sage: K.<a> = NumberField(x^2-x+7)
            sage: a.parts()
            (1/2, 3/2)
            sage: a._coefficients()
            [0, 1]
        """
        cdef Rational ad = <Rational>PY_NEW(Rational), bd = <Rational>PY_NEW(Rational)
        if mpz_cmp_ui(self.a, 0) == 0:
            mpq_set_ui(ad.value, 0, 1)
        else:
            mpz_set(mpq_numref(ad.value), self.a)
            mpz_set(mpq_denref(ad.value), self.denom)
            mpq_canonicalize(ad.value)
        if mpz_cmp_ui(self.b, 0) == 0:
            mpq_set_ui(bd.value, 0, 1)
        else:
            mpz_set(mpq_numref(bd.value), self.b)
            mpz_set(mpq_denref(bd.value), self.denom)
            mpq_canonicalize(bd.value)

        return ad, bd

    cdef bint is_sqrt_disc(self):
        return mpz_cmp_ui(self.denom, 1)==0 and mpz_cmp_ui(self.a,0)==0 and mpz_cmp_ui(self.b, 1)==0


#########################################################
# Arithmatic
#########################################################

    cdef void _reduce_c_(self):
        cdef mpz_t gcd
        mpz_init(gcd)
        mpz_gcd(gcd, self.a, self.denom)
        mpz_gcd(gcd, gcd, self.b)
        if mpz_cmp_si(gcd, 1): # != 0 (i.e. it is not 1)
            mpz_divexact(self.a, self.a, gcd)
            mpz_divexact(self.b, self.b, gcd)
            mpz_divexact(self.denom, self.denom, gcd)
        mpz_clear(gcd)

    cdef ModuleElement _add_c_impl(self, ModuleElement other_m):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2-5)
            sage: K.discriminant()
            5
            sage: a+a
            2*a
            sage: s = (a+2)/6; s
            1/6*a + 1/3
            sage: s+a
            7/6*a + 1/3
            sage: s+10
            1/6*a + 31/3
            sage: s+(2*a+5)/7
            19/42*a + 22/21
            sage: s+(1+a)/2
            2/3*a + 5/6
            sage: s+(1+a)/8
            7/24*a + 11/24
            sage: s+(a+5)/6
            1/3*a + 7/6
            sage: (a/3+2/3) + (2*a/3+1/3)
            a + 1
        """
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        cdef mpz_t gcd, tmp
        res.disc = self.disc
        if mpz_cmp(self.denom, other.denom) == 0:
            mpz_add(res.a, self.a, other.a)
            mpz_add(res.b, self.b, other.b)
            mpz_set(res.denom, self.denom)
        else:
            mpz_init(gcd)
            mpz_gcd(gcd, self.denom, other.denom)
            if mpz_cmp_ui(gcd, 1) == 0:
                mpz_mul(res.a, self.a, other.denom)
                mpz_addmul(res.a, self.denom, other.a)
                mpz_mul(res.b, self.b, other.denom)
                mpz_addmul(res.b, self.denom, other.b)
                mpz_mul(res.denom, self.denom, other.denom)
            else:
                mpz_init(tmp)
                mpz_divexact(tmp, other.denom, gcd)
                mpz_mul(res.a, self.a, tmp)
                mpz_mul(res.b, self.b, tmp)
                mpz_divexact(tmp, self.denom, gcd)
                mpz_addmul(res.a, other.a, tmp)
                mpz_addmul(res.b, other.b, tmp)
                mpz_mul(res.denom, other.denom, tmp)
                mpz_clear(tmp)
            mpz_clear(gcd)
        res._reduce_c_()
        return res

    cdef ModuleElement _sub_c_impl(self, ModuleElement other_m):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2-13)
            sage: b = (a-3)/10; b
            1/10*a - 3/10
            sage: b-1
            1/10*a - 13/10
            sage: b-a
            -9/10*a - 3/10
            sage: b-1/2
            1/10*a - 4/5
            sage: b-a/15
            1/30*a - 3/10
        """
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        cdef mpz_t gcd, tmp
        res.disc = self.disc
        if mpz_cmp(self.denom, other.denom) == 0:
            mpz_sub(res.a, self.a, other.a)
            mpz_sub(res.b, self.b, other.b)
            mpz_set(res.denom, self.denom)
        else:
            mpz_init(gcd)
            mpz_gcd(gcd, self.denom, other.denom)
            if mpz_cmp_ui(gcd, 1) == 0:
                mpz_mul(res.a, self.a, other.denom)
                mpz_submul(res.a, self.denom, other.a)
                mpz_mul(res.b, self.b, other.denom)
                mpz_submul(res.b, self.denom, other.b)
                mpz_mul(res.denom, self.denom, other.denom)
            else:
                mpz_init(tmp)
                mpz_divexact(tmp, other.denom, gcd)
                mpz_mul(res.a, self.a, tmp)
                mpz_mul(res.b, self.b, tmp)
                mpz_divexact(tmp, self.denom, gcd)
                mpz_submul(res.a, other.a, tmp)
                mpz_submul(res.b, other.b, tmp)
                mpz_mul(res.denom, other.denom, tmp)
                mpz_clear(tmp)
            mpz_clear(gcd)
        res._reduce_c_()
        return res

    def __neg__(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+163)
            sage: -a
            -a
            sage: -(a+4)
            -a - 4
            sage: b = (a-3)/2
            sage: -b
            -1/2*a + 3/2
        """
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        res.disc = self.disc
        mpz_neg(res.a, self.a)
        mpz_neg(res.b, self.b)
        mpz_set(res.denom, self.denom)
        return res

    cdef RingElement _mul_c_impl(self, RingElement other_m):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+23)
            sage: a*a
            -23
            sage: (a+1)*(a-1)
            -24
            sage: (a+1)*(a+2)
            3*a - 21
            sage: (a+1)/2 * (a+2)
            3/2*a - 21/2
            sage: (a+1)/2 * (a+2)/3
            1/2*a - 7/2
            sage: (2*a+4) * (3*a)/2
            6*a - 69

        Verify Karatsuba
            sage: K.<a> = NumberField(x^2-41)
            sage: (10^1000 * (a+1)) * K(2+3*a) == 10^1000 * ((a+1) * K(2+3*a))
            True
        """
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        res.disc = self.disc
        cdef mpz_t tmp

        if mpz_size(self.a) + mpz_size(self.b) < 8: # could I use a macro instead?
            # Do it the traditional way
            mpz_mul(res.a, self.b, other.b)
            mpz_mul(res.a, res.a, self.disc.value)
            mpz_addmul(res.a, self.a, other.a)

            mpz_mul(res.b, self.a, other.b)
            mpz_addmul(res.b, self.b, other.a)

        else:
            # Karatsuba
            _sig_on
            mpz_init(tmp)
            mpz_add(res.a, self.a, self.b) # using res.a as tmp
            mpz_add(tmp, other.a, other.b)
            mpz_mul(res.b, res.a, tmp) # res.b = (self.a + self.b)(other.a + other.b)

            mpz_mul(res.a, self.a, other.a)
            mpz_sub(res.b, res.b, res.a)
            mpz_mul(tmp, self.b, other.b)
            mpz_sub(res.b, res.b, tmp)
            mpz_mul(tmp, tmp, self.disc.value)
            mpz_add(res.a, res.a, tmp)
            mpz_clear(tmp)
            _sig_off


        mpz_mul(res.denom, self.denom, other.denom)
        res._reduce_c_()

        return res

    cdef ModuleElement _rmul_c_impl(self, RingElement _c):
        """
        EXAMPLE:
            sage: K.<a> = NumberField(x^2+43)
            sage: (1+a)*3
            3*a + 3
        """
        cdef Rational c = <Rational>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        res.disc = self.disc
        mpz_mul(res.a, self.a, mpq_numref(c.value))
        mpz_mul(res.b, self.b, mpq_numref(c.value))
        mpz_mul(res.denom, self.denom, mpq_denref(c.value))
        res._reduce_c_()
        return res

    cdef ModuleElement _lmul_c_impl(self, RingElement _c):
        """
        EXAMPLE:
            sage: K.<a> = NumberField(x^2+43)
            sage: 5*(a-1/5)
            5*a - 1
        """
        cdef Rational c = <Rational>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        res.disc = self.disc
        mpz_mul(res.a, self.a, mpq_numref(c.value))
        mpz_mul(res.b, self.b, mpq_numref(c.value))
        mpz_mul(res.denom, self.denom, mpq_denref(c.value))
        res._reduce_c_()
        return res

    cdef RingElement _div_c_impl(self, RingElement other):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2-5)
            sage: 2/a
            2/5*a
            sage: (a+2)/(a+1)
            1/4*a + 3/4
            sage: (a+1)*(a+2)/(a+1)
            a + 2
            sage: (a+1/3)*(5*a+2/7)/(a+1/3)
            5*a + 2/7
        """
        return self * ~other

    def __invert__(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2-5)
            sage: ~a
            1/5*a
            sage: ~(a+1)
            1/4*a - 1/4
            sage: (a-1)*(a+1)
            4
            sage: b = ~(5*a-3); b
            5/116*a + 3/116
            sage: b*(5*a-3)
            1
            sage: b = ~((3*a-2)/7); b
            21/41*a + 14/41
            sage: (3*a-2)/7 * b
            1
        """
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        res.disc = self.disc
        cdef mpz_t tmp, gcd
        mpz_init(tmp)
        mpz_init(gcd)

        mpz_gcd(gcd, self.a, self.b)
        if mpz_cmp_si(gcd, 1): # != 0 (i.e. it is not 1)
            # cancel out g (g(a'-b'd)) / (g^2 (a'^2-b'^2d^2))
            mpz_divexact(res.a, self.a, gcd)
            mpz_divexact(res.b, self.b, gcd)
            mpz_neg(res.b, res.b)
        else:
            mpz_set(res.a, self.a)
            mpz_neg(res.b, self.b)

        mpz_pow_ui(res.denom, res.a, 2)
        mpz_pow_ui(tmp, res.b, 2)
        mpz_mul(tmp, tmp, self.disc.value)
        mpz_sub(res.denom, res.denom, tmp)
        # need to multiply the leftover g back in
        mpz_mul(res.denom, res.denom, gcd)

        mpz_mul(res.a, res.a, self.denom)
        mpz_mul(res.b, res.b, self.denom)

        mpz_clear(tmp)
        mpz_clear(gcd)

        res._reduce_c_()
        return res

    cdef NumberFieldElement conjugate_c(self):
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new()
        res.disc = self.disc
        mpz_set(res.a, self.a)
        mpz_neg(res.b, self.b)
        mpz_set(res.denom, self.denom)
        return res

#################################################################################
# We must override everything that makes uses of self.__numerator/__denominator
#################################################################################

    cdef int _cmp_c_impl(self, Element _right) except -2:
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+163)
            sage: K(1/2)==1/2
            True
            sage: a == 1/2
            False
            sage: 2+a == a+2
            True
        """
        cdef NumberFieldElement_quadratic right = _right
        return not mpz_cmp(self.a, right.a)==0  \
            or not mpz_cmp(self.b, right.b)==0  \
            or not mpz_cmp(self.denom, right.denom) == 0

    def __nonzero__(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+163)
            sage: not a
            False
            sage: not (a-a)
            True
        """
        return mpz_cmp_ui(self.a, 0) != 0 or mpz_cmp_ui(self.b, 0) != 0

    def _integer_(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+163)
            sage: (a+1-a)._integer_()
            1
            sage: (a+1/2-a)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce 1/2 to an integer
        """
        cdef Integer res
        if mpz_cmp_ui(self.b, 0) != 0 or mpz_cmp_ui(self.denom, 1) != 0:
            raise TypeError, "Unable to coerce %s to an integer"%self
        else:
            res = PY_NEW(Integer)
            mpz_set(res.value, self.a)
            return res

    def _rational_(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+163)
            sage: (a+1/2-a)._rational_()
            1/2
            sage: (a+1/2)._rational_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce a + 1/2 to a rational
        """
        cdef Rational res
        if mpz_cmp_ui(self.b, 0)!=0:
            raise TypeError, "Unable to coerce %s to a rational"%self
        else:
            res = <Rational>PY_NEW(Rational)
            mpz_set(mpq_numref(res.value), self.a)
            mpz_set(mpq_denref(res.value), self.denom)
            mpq_canonicalize(res.value)
            return res

    def _coefficients(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+41)
            sage: a._coefficients()
            [0, 1]
            sage: K.<a> = NumberField(x^2+x+41)
            sage: a._coefficients()
            [0, 1]
            sage: b = 3*a+1/5
            sage: b._coefficients()
            [1/5, 3]
        """
        # In terms of the generator...
        cdef NumberFieldElement_quadratic gen = self._parent.gen() # should this be cached?
        cdef Rational const = <Rational>PY_NEW(Rational), lin = <Rational>PY_NEW(Rational)
        ad, bd = self.parts()
        if not self:
            return []
        if not bd:
            return [ad]
        if gen.is_sqrt_disc():
            return [ad,bd]
        else:
            alpha, beta = gen.parts()
            scale = bd/beta
            return [ad - scale*alpha, scale]

    def denominator(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+x+41)
            sage: a.denominator()
            1
            sage: b = (2*a+1)/6
            sage: b.denominator()
            6
        """
        # In terms of the generator...
        cdef NumberFieldElement_quadratic gen = self._parent.gen() # should this be cached?
        cdef Integer denom
        if gen.is_sqrt_disc():
            denom = PY_NEW(Integer)
            mpz_set(denom.value, self.denom)
            return denom
        else:
            const, lin = self._coefficients()
            return const.denominator().lcm(lin.denominator())

    cdef bint is_rational_c(self):
        return mpz_cmp_ui(self.b, 0) == 0

#########################################################
# Some things are so much easier to compute
#########################################################

    def trace(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+x+41)
            sage: a.trace()
            -1
            sage: a.matrix()
            [  0   1]
            [-41  -1]

        The trace is additive:
            sage: K.<a> = NumberField(x^2+7)
            sage: (a+1).trace()
            2
            sage: K(3).trace()
            6
            sage: (a+4).trace()
            8
        """
        # trace = 2*self.a / self.denom
        cdef Rational res = <Rational>PY_NEW(Rational)
        if mpz_odd_p(self.denom):
            mpz_mul_2exp(mpq_numref(res.value), self.a, 1)
            mpz_set(mpq_denref(res.value), self.denom)
        else:
            mpz_set(mpq_numref(res.value), self.a)
            mpz_divexact_ui(mpq_denref(res.value), self.denom, 2)
        return res

    def norm(self):
        """
        EXAMPLES:
            sage: K.<a> = NumberField(x^2-x+3)
            sage: a.norm()
            3
            sage: a.matrix()
            [ 0  1]
            [-3  1]
            sage: K.<a> = NumberField(x^2+5)
            sage: (1+a).norm()
            6

        The norm is multiplicative:
            sage: K.<a> = NumberField(x^2-3)
            sage: a.norm()
            -3
            sage: K(3).norm()
            9
            sage: (3*a).norm()
            -27
        """
        # norm = (a^2 - d b^2) / self.denom^2
        cdef Rational res = <Rational>PY_NEW(Rational)
        mpz_pow_ui(mpq_numref(res.value), self.a, 2)
        mpz_pow_ui(mpq_denref(res.value), self.b, 2) # use as temp
        mpz_mul(mpq_denref(res.value), mpq_denref(res.value), self.disc.value)
        mpz_sub(mpq_numref(res.value), mpq_numref(res.value), mpq_denref(res.value))
        mpz_pow_ui(mpq_denref(res.value), self.denom, 2)
        mpq_canonicalize(res.value)
        return res

    def charpoly(self, var='x'):
        r"""
        The characteristic polynomial of this element over $\Q$.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2-x+13)
            sage: a.charpoly()
            x^2 - x + 13
            sage: b = 3-a/2
            sage: f = b.charpoly(); f
            x^2 - 11/2*x + 43/4
            sage: f(b)
            0
        """
        R = QQ[var]
        return R([self.norm(), -self.trace(), 1])

    def minopily(self, var='x'):
        r"""
        The minimal polynomial of this element over $\Q$.
        EXAMPLES:
            sage: K.<a> = NumberField(x^2+13)
            sage: a.minpoly()
            x^2 + 13
            sage: (a+1/2-a).minpoly()
            x - 1/2
        """

        if self.is_rational_c():
            R = QQ[var]
            return R([-self._rational_(), 1])
        else:
            return self.charpoly()
