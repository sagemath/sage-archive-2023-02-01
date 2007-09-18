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

AUTHORS:
    -- Robert Bradshaw (2007-09): Initial version
"""

include '../../ext/interrupt.pxi'
include "../../ext/stdsage.pxi"
cdef extern from *:
    # TODO: move to stdsage.pxi
    object PY_NEW_SAME_TYPE(object o)
    bint PY_TYPE_CHECK_EXACT(object o, object type)

cdef object QQ, ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ



cdef class NumberFieldElement_quadratic(NumberFieldElement_absolute):
    # (a + b sqrt(disc)) / denom
    def __init__(self, parent, f):
        NumberFieldElement_absolute.__init__(self, parent, f) # this is heavy, can I avoid it for f in Rational, tuple?
        self.disc = parent.discriminant()
        cdef Integer a, b, denom
        cdef Rational ad, bd

        cdef Py_ssize_t fdeg
        cdef mpz_t tmp
        cdef NumberFieldElement_quadratic gen

        if PY_TYPE_CHECK_EXACT(f, Rational):
            mpz_set(self.a, mpq_numref((<Rational>f).value))
            mpz_set_ui(self.b, 0)
            mpz_set(self.denom, mpq_denref((<Rational>f).value))

        elif PY_TYPE_CHECK_EXACT(f, tuple) and len(f) == 2:
            ad, bd = f
            mpz_set(self.a, a.value)
            mpz_set(self.b, b.value)
            mpz_lcm(self.denom, mpq_denref(ad.value), mpq_denref(bd.value))
            mpz_divexact(self.a, self.denom, mpq_denref(ad.value))
            mpz_mul(self.a, self.a, mpq_numref(ad.value))
            mpz_divexact(self.b, self.denom, mpq_denref(bd.value))
            mpz_mul(self.b, self.b, mpq_numref(bd.value))

        elif PY_TYPE_CHECK_EXACT(f, tuple) and len(f) == 3:
            a, b, denom = f
            mpz_set(self.a, a.value)
            mpz_set(self.b, b.value)
            mpz_set(self.denom, denom.value)
            self._reduce_c()

        else:
            # poly in gen (which may not be sqrt(d))
            self._ntl_coeff_as_mpz(&self.a, 0)
            if mpz_cmp_ui(self.a, 0):
                self._ntl_coeff_as_mpz(&self.b, 1)
                if mpz_cmp_ui(self.b, 0):
                    gen = parent.gen() # should this be cached?
                    mpz_init(tmp)
                    mpz_mul(tmp, self.b, gen.a)
                    mpz_clear(tmp)
                    mpz_add(self.a, self.a, tmp)
                    mpz_mul(self.b, self.b, gen.b)
                self._ntl_denom_as_mpz(&self.denom)
            else:
                mpz_set_ui(self.b, 0)
                mpz_set_ui(self.denom, 1)

    def __new__(self):
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.denom)

    def __dealloc__(self):
        mpz_clear(self.a)
        mpz_clear(self.b)
        mpz_clear(self.denom)

    def parts(self):
        cdef Rational ad = <Rational>PY_NEW(Rational), bd = <Rational>PY_NEW(Rational)
        mpz_set(mpq_numref(ad.value), self.a)
        mpz_set(mpq_denref(ad.value), self.denom)
        mpz_set(mpq_numref(bd.value), self.b)
        mpz_set(mpq_denref(bd.value), self.denom)
        return ad, bd

    cdef bint is_sqrt_disc(self):
        return mpz_cmp_ui(self.denom, 1)==0 and mpz_cmp_ui(self.a,0)==0 and mpz_cmp_ui(self.b, 1)==0


#########################################################
# Arithmatic
#########################################################

    cdef void _reduce_c_(self):
        cdef mpz_t gcd
        mpz_init(gcd)
        mpz_gcd(gcd, self.a, self.b)
        mpz_gcd(gcd, gcd, self.denom)
        if mpz_cmp_si(gcd, 1): # != 0 (i.e. it is not 1)
            mpz_divexact(self.a, self.a, gcd)
            mpz_divexact(self.b, self.b, gcd)
            mpz_divexact(self.denom, self.denom, gcd)
        mpz_clear(gcd)

    cdef ModuleElement _add_c_impl(self, ModuleElement other_m):
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
        res.disc = self.disc
        mpz_add(res.a, self.a, other.a)
        mpz_add(res.b, self.b, other.b)
        mpz_set(res.denom, self.denom)
        res._reduce_c_()
        return res

    cdef ModuleElement _sub_c_impl(self, ModuleElement other_m):
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
        res.disc = self.disc
        mpz_sub(res.a, self.a, other.a)
        mpz_sub(res.b, self.b, other.b)
        mpz_set(res.denom, self.denom)
        res._reduce_c_()
        return res

    def __neg__(self):
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
        res.disc = self.disc
        mpz_neg(res.a, self.a)
        mpz_neg(res.b, self.b)
        mpz_set(res.denom, self.denom)
        return res

    cdef RingElement _mul_c_impl(self, RingElement other_m):
        cdef NumberFieldElement_quadratic other = <NumberFieldElement_quadratic>other_m
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
        res.disc = self.disc
        cdef mpz_t tmp
        mpz_init(tmp)

        if mpz_size(self.a) + mpz_size(self.b) < 8: # could I use a macro instead?
            # Do it the traditional way
            mpz_mul(res.a, self.a, other.a)
            mpz_mul(tmp, self.b, other.b)
            mpz_mul(tmp, tmp, self.disc.value)
            mpz_add(res.a, res.a, tmp)

            mpz_mul(res.b, self.a, other.b)
            mpz_mul(tmp, self.a, other.b)
            mpz_add(res.b, res.b, tmp)

        else:
            # Karatsuba
            _sig_on
            mpz_add(res.a, self.a, other.a) # using res.a as tmp
            mpz_add(tmp, self.b, other.b)
            mpz_mul(res.b, res.a, tmp) # res.b = (self.a + other.a)(self.b + other.b)

            mpz_mul(res.a, self.a, other.a)
            mpz_sub(res.b, res.b, res.a)
            mpz_mul(tmp, self.b, other.b)
            mpz_sub(res.b, res.b, tmp)
            mpz_mul(tmp, tmp, self.disc.value)
            mpz_add(res.a, tmp, tmp)
            _sig_off

        mpz_clear(tmp)

        mpz_mul(res.denom, self.denom, other.denom)
        res._reduce_c_()

        return res

    cdef ModuleElement _rmul_c_impl(self, RingElement _c):
        cdef Rational c = <Rational>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
        res.disc = self.disc
        mpz_mul(res.a, self.a, mpq_numref(c.value))
        mpz_mul(res.b, self.b, mpq_numref(c.value))
        mpz_mul(res.denom, self.denom, mpq_denref(c.value))
        res._reduce_c_()
        return res

    cdef ModuleElement _lmul_c_impl(self, RingElement _c):
        cdef Rational c = <Rational>_c
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
        res.disc = self.disc
        mpz_mul(res.a, self.a, mpq_numref(c.value))
        mpz_mul(res.b, self.b, mpq_numref(c.value))
        mpz_mul(res.denom, self.denom, mpq_denref(c.value))
        res._reduce_c_()
        return res

    cdef RingElement _div_c_impl(self, RingElement other):
        return self * ~other

    def __invert__(self):
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
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

        mpz_mul(res.denom, res.denom, self.denom)

        mpz_clear(tmp)
        mpz_clear(gcd)

        res._reduce_c_()
        return res

    cdef NumberFieldElement conjugate_c(self):
        cdef NumberFieldElement_quadratic res = <NumberFieldElement_quadratic>self._new_c()
        res.disc = self.disc
        mpz_set(res.a, self.a)
        mpz_neg(res.b, self.b)
        mpz_set(res.denom, self.denom)
        return res

#################################################################################
# We must override everything that makes uses of self.__numerator/__denominator
#################################################################################

    cdef int _cmp_c_impl(self, Element _right) except -2:
        cdef NumberFieldElement_quadratic right = _right
        return not mpz_cmp(self.a, right.a)==0  \
            or not mpz_cmp(self.b, right.b)==0  \
            or not mpz_cmp(self.denom, right.denom) == 0

    def __nonzero__(self):
        return mpz_cmp_ui(self.a, 0) != 0 or mpz_cmp_ui(self.b, 0) != 0

    def _integer_(self):
        cdef Integer res
        if mpz_cmp_ui(self.b, 0) != 0 or mpz_cmp_ui(self.denom, 1) != 0:
            raise TypeError, "Unable to coerce %s to an integer"%self
        else:
            res = PY_NEW(Integer)
            mpz_set(res.value, self.a)
            return res

    def _rational_(self):
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
        # In terms of the generator...
        cdef NumberFieldElement_quadratic gen = self._parent.gen() # should this be cached?
        cdef Rational const = <Rational>PY_NEW(Rational), lin = <Rational>PY_NEW(Rational)
        if gen.is_sqrt_disc():
            # gen = sqrt(disc)
            ad, bd = self.parts()
            return [ad,bd]
        else:
            alpha, beta = gen.parts()
            scale = bd/beta
            return [ad - scale*alpha, scale]

    def denominator(self):
        cdef Integer denom = PY_NEW(Integer)
        mpz_set(denom.value, self.denom)
        return denom

    cdef bint is_rational_c(self):
        return mpz_cmp_ui(self.b, 0) == 0

#########################################################
# Some things are so much easier to compute
#########################################################

    def trace(self):
        # trace = 2*a
        cdef Rational res = <Rational>PY_NEW(Rational)
        if mpz_odd_p(self.denom):
            mpz_mul_2exp(mpq_numref(res.value), self.a, 1)
            mpz_set(mpq_denref(res.value), self.denom)
        else:
            mpz_set(mpq_numref(res.value), self.a)
            mpz_divexact_ui(mpq_denref(res.value), self.denom, 2)

    def norm(self):
        # norm = a^2 - d b^2
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

        We compute the charpoly of cube root of $2$.

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2-7)
            sage: a.charpoly('x')
            x^2 - 7

        """
        R = QQ[var]
        return R([self.norm(), -self.trace(), 1])
