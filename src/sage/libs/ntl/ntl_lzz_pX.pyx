"""
ntl_lzz_pX.pyx

Wraps NTL's zz_pX type for SAGE

AUTHORS:
   - Craig Citro
"""

#*****************************************************************************
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

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
include 'misc.pxi'
include 'decl.pxi'

from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer cimport Integer
from sage.rings.integer_ring cimport IntegerRing_class

from sage.rings.integer_mod cimport IntegerMod_gmp, IntegerMod_int, IntegerMod_int64

from sage.libs.ntl.ntl_lzz_pContext import ntl_zz_pContext
from sage.libs.ntl.ntl_lzz_pContext cimport ntl_zz_pContext_class

ZZ_sage = IntegerRing()

##############################################################################
#
# zz_pX  -- polynomials over the integers modulo p, p small
#
##############################################################################

cdef class ntl_zz_pX:
    r"""
    The class \class{zz_pX} implements polynomial arithmetic modulo $p$,
    for p smaller than a machine word.

    Polynomial arithmetic is implemented using the FFT, combined with
    the Chinese Remainder Theorem.  A more detailed description of the
    techniques used here can be found in [Shoup, J. Symbolic
    Comp. 20:363-397, 1995].

    Small degree polynomials are multiplied either with classical
    or Karatsuba algorithms.
    """
    # See ntl_zz_pX.pxd for definition of data members
    def __init__(self, ls=[], modulus=None):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([1,2,5,-9])
            sage: f
            [1 2 5 11]
            sage: g = ntl.zz_pX([0,0,0]); g
            []
            sage: g[10]=5
            sage: g
            [0 0 0 0 0 0 0 0 0 0 5]
            sage: g[10]
            5
            sage: f = ntl.zz_pX([10^30+1, 10^50+1], 100); f
            [1 1]
        """
        cdef long n
        cdef Py_ssize_t i
        cdef long temp

        if PY_TYPE_CHECK( modulus, ntl_zz_pContext_class ):
            self.c = <ntl_zz_pContext_class>modulus
            p_sage = Integer(self.c.p)
        elif PY_TYPE_CHECK( modulus, Integer ):
            if modulus > NTL_SP_BOUND:
                raise ValueError, "Modulus (=%s) is too big" % modulus
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(mpz_get_si((<Integer>modulus).value))
            p_sage = modulus
        elif PY_TYPE_CHECK( modulus, long ):
            if modulus > NTL_SP_BOUND:
                raise ValueError, "Modulus (=%s) is too big" % modulus
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(modulus)
            p_sage = Integer(self.c.p)
        elif modulus is None:
            raise ValueError, "You must specify a modulus."
        else:
            try:
                modulus = long(modulus)
                if modulus > NTL_SP_BOUND:
                    raise ValueError, "Modulus (=%s) is too big" % modulus
            except:
                raise ValueError, "%s (type %s) is not a valid modulus." % (modulus, type(modulus))
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(modulus)
            p_sage = Integer(self.c.p)

        ## now that we've determined the modulus, set that modulus.
        self.c.restore_c()

        n = len(ls)
        if (n == 0):
            ## the 0 polynomial is just the empty list;
            ## so in this case, we're done.
            return

        self.x.SetMaxLength(n+1)

        for i from 0 <= i < n:
            a = ls[i]

            if PY_TYPE_CHECK(a, IntegerMod_int):
                if (self.c.p == (<IntegerMod_int>a).__modulus.int32): ## this is slow
                    zz_pX_SetCoeff_long(self.x, i, (<IntegerMod_int>a).ivalue)
                else:
                    raise ValueError, \
                          "Mismatched modulus for converting to zz_pX."
            elif PY_TYPE_CHECK(a, IntegerMod_int64):
                if (self.c.p == (<IntegerMod_int64>a).__modulus.int64): ## this is slow
                    zz_pX_SetCoeff_long(self.x, i, (<IntegerMod_int64>a).ivalue)
                else:
                    raise ValueError, \
                          "Mismatched modulus for converting to zz_pX."
            elif PY_TYPE_CHECK(a, IntegerMod_gmp):
                if (p_sage == (<IntegerMod_gmp>a).__modulus.sageInteger): ## this is slow
                    zz_pX_SetCoeff_long(self.x, i, mpz_get_si((<IntegerMod_gmp>a).value))
                else:
                    raise ValueError, \
                          "Mismatched modulus for converting to zz_pX."
            elif PY_TYPE_CHECK(a, Integer):
                zz_pX_SetCoeff_long(self.x, i, mpz_fdiv_ui((<Integer>a).value, self.c.p))
            elif PY_TYPE_CHECK(a, int):
                ## we're lucky that python int is no larger than long
                temp = a
                zz_pX_SetCoeff_long(self.x, i, temp%self.c.p)
            else:
                a = Integer(a)
                zz_pX_SetCoeff_long(self.x, i, mpz_fdiv_ui((<Integer>a).value, self.c.p))

##        zz_pX_SetCoeff_long(self.x, n, <long>1)

        return

    def __reduce__(self):
        return make_zz_pX, (self.list(), self.c)

    def __repr__(self):
        return str(self.list())

    cdef ntl_zz_pX _new(self):
        cdef ntl_zz_pX y
        y = PY_NEW(ntl_zz_pX)
        y.c = self.c
        return y

    def __add__(ntl_zz_pX self, other):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: ntl.zz_pX(range(5)) + ntl.zz_pX(range(6))
            [0 2 4 6 8 5]
        """
        if not PY_TYPE_CHECK(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._add_c_impl(other)

    cdef ntl_zz_pX _add_c_impl(ntl_zz_pX self, ntl_zz_pX other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        zz_pX_add(y.x, self.x, other.x)
        _sig_off
        return y

    def __sub__(ntl_zz_pX self, other):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: ntl.zz_pX(range(5)) - ntl.zz_pX(range(6))
            [0 0 0 0 0 15]
        """
        if not PY_TYPE_CHECK(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._sub_c_impl(other)

    cdef ntl_zz_pX _sub_c_impl(ntl_zz_pX self, ntl_zz_pX other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        zz_pX_sub(y.x, self.x, other.x)
        _sig_off
        return y

    def __mul__(ntl_zz_pX self, other):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: ntl.zz_pX(range(5)) * ntl.zz_pX(range(6))
            [0 0 1 4 10 0 10 14 11]
        """
        if not PY_TYPE_CHECK(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._mul_c_impl(other)

    cdef ntl_zz_pX _mul_c_impl(ntl_zz_pX self, ntl_zz_pX other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        zz_pX_mul(y.x, self.x, other.x)
        _sig_off
        return y

    def __div__(ntl_zz_pX self, other):
        """
        Compute quotient self / other, if the quotient is a polynomial.
        Otherwise an Exception is raised.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.zz_pX([1,2,3]) * ntl.zz_pX([4,5])**2
            sage: g = ntl.zz_pX([4,5])
            sage: f/g
            [4 13 5 15]
            sage: ntl.zz_pX([1,2,3]) * ntl.zz_pX([4,5])
            [4 13 5 15]

            sage: f = ntl.zz_pX(range(10)); g = ntl.zz_pX([-1,0,1])
            sage: f/g
            Traceback (most recent call last):
            ...
            ArithmeticError: self (=[0 1 2 3 4 5 6 7 8 9]) is not divisible by other (=[16 0 1])
        """
        if not PY_TYPE_CHECK(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._div_c_impl(other)

    cdef ntl_zz_pX _div_c_impl(ntl_zz_pX self, ntl_zz_pX other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        _sig_on
        cdef long divisible
        cdef ntl_zz_pX q
        q = self._new()
        divisible = zz_pX_divide(q.x, self.x, other.x)
        if not divisible:
            raise ArithmeticError, "self (=%s) is not divisible by other (=%s)"%(self, other)
        ## FIXME
        return q

    def __mod__(ntl_zz_pX self, other):
        """
        Given polynomials a, b in ZZ[X], there exist polynomials q, r
        in QQ[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns q if q lies in ZZ[X], and otherwise raises an
        Exception.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.zz_pX([2,4,6]); g = ntl.zz_pX([2])
            sage: f % g   # 0
            []

            sage: f = ntl.zz_pX(range(10)); g = ntl.zz_pX([-1,0,1])
            sage: f % g
            [3 8]
        """
        if not PY_TYPE_CHECK(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._mod_c_impl(other)

    cdef ntl_zz_pX _mod_c_impl(ntl_zz_pX self, ntl_zz_pX other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        zz_pX_mod(y.x, self.x, other.x)
        _sig_off
        return y

    def __pow__(ntl_zz_pX self, long n, ignored):
        """
        Return the n-th nonnegative power of self.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: g = ntl.zz_pX([-1,0,1])
            sage: g**10
            [1 0 10 0 5 0 0 0 10 0 8 0 10 0 0 0 5 0 10 0 1]
        """
        if n < 0:
            raise ValueError, "Only positive exponents allowed."
        _sig_on
        cdef ntl_zz_pX y = self._new()
        zz_pX_power(y.x, self.x, n)
        _sig_off
        return y

    def quo_rem(ntl_zz_pX self, ntl_zz_pX right):
        """
        Returns the quotient and remander when self is devided by right.

        Specifically, this return r, q such that $self = q * right + r$

        EXAMPLES:
            sage: f = ntl.zz_pX(range(7), 19)
            sage: g = ntl.zz_pX([2,4,6], 19)
            sage: f // g
            [1, 1, 15, 16, 1]
            sage: f % g
            [17, 14]
            sage: f.quo_rem(g)
            ([1, 1, 15, 16, 1], [17, 14])
            sage: (f // g) * g + f % g
            [0, 1, 2, 3, 4, 5, 6]
        """
        cdef ntl_zz_pX q = self._new()
        cdef ntl_zz_pX r = self._new()
        _sig_on
        self.c.restore_c()
        zz_pX_divrem(q.x, r.x, self.x, right.x)
        _sig_off
        return q, r

    def __floordiv__(ntl_zz_pX self, ntl_zz_pX right):
        """
        Returns the whole part of $self / right$.

        EXAMPLE:
            sage: f = ntl.zz_pX(range(10), 19); g = ntl.zz_pX([1]*5, 19)
            sage: f // g
            [8, 18, 18, 18, 18, 9]

        """
        _sig_on
        self.c.restore_c()
        cdef ntl_zz_pX q = self._new()
        zz_pX_div(q.x, self.x, right.x)
        _sig_off
        return q

    def __lshift__(ntl_zz_pX self, long n):
        """
        Shifts this polynomial to the left, which is multiplication by $x^n$.

        EXAMPLE:
            sage: f = ntl.zz_pX([2,4,6], 17)
            sage: f << 2
            [0, 0, 2, 4, 6]
        """
        cdef ntl_zz_pX r = self._new()
        zz_pX_lshift(r.x, self.x, n)
        return r

    def __rshift__(ntl_zz_pX self, long n):
        """
        Shifts this polynomial to the right, which is division by $x^n$ (and truncation).

        EXAMPLE:
            sage: f = ntl.zz_pX([1,2,3], 17)
            sage: f >> 2
            [3]
        """
        cdef ntl_zz_pX r = self._new()
        zz_pX_rshift(r.x, self.x, n)
        return r

    def diff(self):
        """
        The formal derivative of self.

        EXAMPLE:
            sage: f = ntl.zz_pX(range(10), 17)
            sage: f.diff()
            [1, 4, 9, 16, 8, 2, 15, 13, 13]
        """
        cdef ntl_zz_pX r = self._new()
        zz_pX_diff(r.x, self.x)
        return r

    def reverse(self):
        """
        Returns self with coefficients reversed, i.e. $x^n self(x^{-n})$.

        EXAMPLE:
            sage: f = ntl.zz_pX([2,4,6], 17)
            sage: f.reverse()
            [6, 4, 2]
        """
        cdef ntl_zz_pX r = self._new()
        zz_pX_reverse(r.x, self.x)
        return r

    def __neg__(self):
        """
        Return the negative of self.
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([2,0,0,1])
            sage: -f
            [18 0 0 19]
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        zz_pX_neg(y.x, self.x)
        _sig_off
        return y

    def __cmp__(ntl_zz_pX self, ntl_zz_pX other):
        """
        Decide whether or not self and other are equal.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([1,2,3])
            sage: g = ntl.zz_pX([1,2,3,0])
            sage: f == g
            True
            sage: g = ntl.zz_pX([0,1,2,3])
            sage: f == g
            False
        """
        if (NTL_zz_pX_DOUBLE_EQUALS(self.x, other.x)):
            return 0
        else:
            return -1

    def list(self):
        """
        Return list of entries as a list of python ints.

        EXAMPLES:
            sage: f = ntl.zz_pX([23, 5,0,1], 10)
            sage: f.list()
            [3, 5, 0, 1]
            sage: type(f.list()[0])
            <type 'int'>
        """
        cdef long i
        return [ zz_p_rep(zz_pX_GetCoeff(self.x, i)) for i from 0 <= i <= zz_pX_deg(self.x) ]

    def degree(self):
        """
        Return the degree of this polynomial.  The degree of the 0
        polynomial is -1.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([5,0,1])
            sage: f.degree()
            2
            sage: f = ntl.zz_pX(range(100))
            sage: f.degree()
            99
            sage: f = ntl.zz_pX()
            sage: f.degree()
            -1
            sage: f = ntl.zz_pX([1])
            sage: f.degree()
            0
        """
        return zz_pX_deg(self.x)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([3,6,9])
            sage: f.leading_coefficient()
            9
            sage: f = ntl.zz_pX()
            sage: f.leading_coefficient()
            0
        """
        return zz_p_rep(zz_pX_LeadCoeff(self.x))

    def constant_term(self):
        """
        Return the constant coefficient of this polynomial.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([3,6,9])
            sage: f.constant_term()
            3
            sage: f = ntl.zz_pX()
            sage: f.constant_term()
            0
        """
        return zz_p_rep(zz_pX_ConstTerm(self.x))

    def square(self):
        """
        Return f*f.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([-1,0,1])
            sage: f*f
            [1 0 15 0 1]
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        zz_pX_sqr(y.x, self.x)
        _sig_off
        return y

    def truncate(self, long m):
        """
        Return the truncation of this polynomial obtained by
        removing all terms of degree >= m.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([1,2,3,4,5])
            sage: f.truncate(3)
            [1 2 3]
            sage: f.truncate(8)
            [1 2 3 4 5]
            sage: f.truncate(1)
            [1]
            sage: f.truncate(0)
            []
            sage: f.truncate(-1)
            []
            sage: f.truncate(-5)
            []
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            zz_pX_trunc(y.x, self.x, m)
        _sig_off
        return y

    def multiply_and_truncate(self, ntl_zz_pX other, long m):
        """
        Return self*other but with terms of degree >= m removed.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([1,2,3,4,5])
            sage: g = ntl.zz_pX([10])
            sage: f.multiply_and_truncate(g, 2)
            [10]
            sage: g.multiply_and_truncate(f, 2)
            [10]
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            zz_pX_MulTrunc(y.x, self.x, other.x, m)
        _sig_off
        return y

    def square_and_truncate(self, long m):
        """
        Return self*self but with terms of degree >= m removed.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([1,2,3,4,5])
            sage: f.square_and_truncate(4)
            [1 4 10]
            sage: (f*f).truncate(4)
            [1 4 10]
        """
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            zz_pX_SqrTrunc(y.x, self.x, m)
        _sig_off
        return y

    def invert_and_truncate(self, long m):
        """
        Compute and return the inverse of self modulo $x^m$.
        The constant term of self must be 1 or -1.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([1,2,3,4,5,6,7])
            sage: f.invert_and_truncate(20)
            [1 18 1 0 0 0 0 8 17 2 13 0 0 0 4 0 17 10 9]
            sage: g = f.invert_and_truncate(20)
            sage: g * f
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 4 4 3]
        """
        if m < 0:
            raise ArithmeticError, "m (=%s) must be positive"%m
        n = self.constant_term()
        if n != 1 and n != -1:
            raise ArithmeticError, \
                  "The constant term of self must be 1 or -1."
        _sig_on
        cdef ntl_zz_pX y
        y = self._new()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            zz_pX_InvTrunc(y.x, self.x, m)
        _sig_off
        return y


    def is_zero(self):
        """
        Return True exactly if this polynomial is 0.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([0,0,0,20])
            sage: f.is_zero()
            True
            sage: f = ntl.zz_pX([0,0,1])
            sage: f
            [0 0 1]
            sage: f.is_zero()
            False
        """
        return zz_pX_IsZero(self.x)

    def is_one(self):
        """
        Return True exactly if this polynomial is 1.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([1,1])
            sage: f.is_one()
            False
            sage: f = ntl.zz_pX([1])
            sage: f.is_one()
            True
        """
        return zz_pX_IsOne(self.x)

    def is_monic(self):
        """
        Return True exactly if this polynomial is monic.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX([2,0,0,1])
            sage: f.is_monic()
            True
            sage: g = f.reverse()
            sage: g.is_monic()
            False
            sage: g
            [1 0 0 2]
            sage: f = ntl.zz_pX([1,2,0,3,0,2])
            sage: f.is_monic()
            False
        """
        # The following line is what we should have.  However, strangely this is *broken*
        # on PowerPC Intel in NTL, so we program around
        # the problem.  (William Stein)
        #return bool(ZZ_pX_is_monic(self.x))

        if zz_pX_IsZero(self.x):
             return False
        return ( zz_p_rep(zz_pX_LeadCoeff(self.x)) == 1 )

    def set_x(self):
        """
        Set this polynomial to the monomial "x".

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX()
            sage: f.set_x()
            sage: f
            [0 1]
            sage: g = ntl.zz_pX([0,1])
            sage: f == g
            True

        Though f and g are equal, they are not the same objects in memory:
            sage: f is g
            False

        """
        zz_pX_SetX(self.x)

    def is_x(self):
        """
        True if this is the polynomial "x".

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.zz_pX()
            sage: f.set_x()
            sage: f.is_x()
            True
            sage: f = ntl.zz_pX([0,1])
            sage: f.is_x()
            True
            sage: f = ntl.zz_pX([1])
            sage: f.is_x()
            False
        """
        return zz_pX_IsX(self.x)

    def clear(self):
        """
        Reset this polynomial to 0.  Changes this polynomial in place.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.zz_pX([1,2,3])
            sage: f
            [1 2 3]
            sage: f.clear()
            sage: f
            []
        """
        zz_pX_clear(self.x)

    def preallocate_space(self, long n):
        """
        Pre-allocate spaces for n coefficients.  The polynomial that f
        represents is unchanged.  This is useful if you know you'll be
        setting coefficients up to n, so memory isn't re-allocated as
        the polynomial grows.  (You might save a millisecond with this
        function.)

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.zz_pX([1,2,3])
            sage: f.preallocate_space(20)
            sage: f
            [1 2 3]
            sage: f[10]=5  # no new memory is allocated
            sage: f
            [1 2 3 0 0 0 0 0 0 0 5]
        """
        _sig_on
        self.x.SetMaxLength(n)
        _sig_off
        return


def make_zz_pX(L, context):
    """
    For unpickling.

    TEST:
        sage: f = ntl.zz_pX(range(16), 12)
        sage: loads(dumps(f)) == f
        True
    """
    return ntl_zz_pX(L, context)
