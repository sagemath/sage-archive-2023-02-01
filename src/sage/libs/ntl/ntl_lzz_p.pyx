"""
ntl_lzz_p.pyx

Wraps NTL's zz_p type for SAGE

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

cdef class ntl_zz_p:
    r"""
    The class \class{zz_p} implements polynomial arithmetic modulo $p$,
    for p smaller than a machine word.

    Polynomial arithmetic is implemented using the FFT, combined with
    the Chinese Remainder Theorem.  A more detailed description of the
    techniques used here can be found in [Shoup, J. Symbolic
    Comp. 20:363-397, 1995].

    Small degree polynomials are multiplied either with classical
    or Karatsuba algorithms.
    """
    # See ntl_zz_p.pxd for definition of data members
    def __init__(self, a=0, modulus=None):
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
        """
        cdef long temp

        if PY_TYPE_CHECK( modulus, ntl_zz_pContext_class ):
            self.c = <ntl_zz_pContext_class>modulus
            p_sage = Integer(self.c.p)
        elif PY_TYPE_CHECK( modulus, Integer ):
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(mpz_get_si((<Integer>modulus).value))
            p_sage = modulus
        elif PY_TYPE_CHECK( modulus, long ):
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(modulus)
            p_sage = Integer(self.c.p)
        elif modulus is None:
            raise ValueError, "You must specify a modulus."
        else:
            try:
                modulus = long(modulus)
            except:
                raise ValueError, "%s is not a valid modulus."%modulus
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(modulus)
            p_sage = Integer(self.c.p)

        ## now that we've determined the modulus, set that modulus.
        self.c.restore_c()

        ## old code to set modulus
##        p_long = <long>p
##        if (p_long > NTL_SP_BOUND):
##            raise ValueError, "Modulus (=%s) is too big"%p
##        zz_p_set_modulus(p_long)

##        p = self.c.p
##        p_sage = Integer(self.c.p)

        if PY_TYPE_CHECK(a, IntegerMod_int):
            if (self.c.p == (<IntegerMod_int>a).__modulus.int32): ## this is slow
                zz_p_set_from_long(self.x, (<IntegerMod_int>a).ivalue)
            else:
                raise ValueError, \
                      "Mismatched modulus for converting to zz_p."

        elif PY_TYPE_CHECK(a, IntegerMod_int64):
            if (self.c.p == (<IntegerMod_int64>a).__modulus.int64): ## this is slow
                zz_p_set_from_long(self.x, (<IntegerMod_int64>a).ivalue)
            else:
                raise ValueError, \
                      "Mismatched modulus for converting to zz_p."

        elif PY_TYPE_CHECK(a, IntegerMod_gmp):
            if (p_sage == (<IntegerMod_gmp>a).__modulus.sageInteger): ## this is slow
                zz_p_set_from_long(self.x, mpz_get_si((<IntegerMod_gmp>a).value))
            else:
                raise ValueError, \
                      "Mismatched modulus for converting to zz_p."

        elif PY_TYPE_CHECK(a, Integer):
            zz_p_set_from_long(self.x, mpz_get_si((<Integer>a).value)%self.c.p)

        elif PY_TYPE_CHECK(a, int):
            ## we're lucky that python int is no larger than long
            temp = a
            zz_p_set_from_long(self.x, temp%self.c.p)
        else:
            a = Integer(a)
            zz_p_set_from_long(self.x, mpz_get_si((<Integer>a).value)%self.c.p)

        return

##    TODO: these are in ntl_ZZ_p.pyx; should I have versions
##    here too?
##    def __new__(self, a=None, modulus=None):
##        zz_p_construct(&self.x)

##    def __dealloc__(self):
##        if <object>self.c is not None:
##            self.c.restore_c()
##        zz_p_destruct(&self.x)

    cdef ntl_zz_p _new(self):
        cdef ntl_zz_p y
        y = PY_NEW(ntl_zz_p)
        y.c = self.c
        return y

    def __reduce__(self):
        raise NotImplementedError

    def __repr__(self):
        return Integer(zz_p_rep(self.x)).__repr__()
##        return ZZ_sage(self.x)._repr_()
##        return str(self.list())

    def __add__(ntl_zz_p self, other):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: ntl.zz_pX(range(5)) + ntl.zz_pX(range(6))
            [0 2 4 6 8 5]
        """
        if not PY_TYPE_CHECK(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._add_c_impl(other)

    cdef ntl_zz_p _add_c_impl(ntl_zz_p self, ntl_zz_p other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        cdef ntl_zz_p y
        y = self._new()
        _sig_on
        zz_p_add(y.x, self.x, other.x)
        _sig_off
        return y

    def __sub__(ntl_zz_p self, other):
        if not PY_TYPE_CHECK(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._sub_c_impl(other)

    cdef ntl_zz_p _sub_c_impl(ntl_zz_p self, ntl_zz_p other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        _sig_on
        cdef ntl_zz_p y
        y = self._new()
        zz_p_sub(y.x, self.x, other.x)
        _sig_off
        return y

    def __mul__(ntl_zz_p self, other):
        if not PY_TYPE_CHECK(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
#        _sig_on
        self.c.restore_c()
#        _sig_off
        return self._mul_c_impl(other)

    cdef ntl_zz_p _mul_c_impl(ntl_zz_p self, ntl_zz_p other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
#        _sig_on
        cdef ntl_zz_p y
        y = self._new()
        zz_p_mul(y.x, self.x, other.x)
#        _sig_off
        return y

    def __div__(ntl_zz_p self, other):
        if not PY_TYPE_CHECK(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError, "arithmetic operands must have the same modulus."
        _sig_on
        self.c.restore_c()
        _sig_off
        return self._div_c_impl(other)

    cdef ntl_zz_p _div_c_impl(ntl_zz_p self, ntl_zz_p other):
        """
        Quick and dirty -- no typechecking, assumes
        that the context is correct, and that the operands have the
        same modulus.
        """
        _sig_on
        cdef ntl_zz_p q
        q = self._new()
        zz_p_div(q.x, self.x, other.x)
        return q

    def __pow__(ntl_zz_p self, long n, ignored):
        """
        Return the n-th nonnegative power of self.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: g = ntl.zz_pX([-1,0,1])
            sage: g**10
            [1 0 10 0 5 0 0 0 10 0 8 0 10 0 0 0 5 0 10 0 1]
        """
        ## FIXME
        if n < 0:
            raise NotImplementedError
        import sage.rings.arith
        return sage.rings.arith.generic_power(self, n, ntl_zz_p([1]))

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
        cdef ntl_zz_p y
        y = self._new()
        self.c.restore_c()
        zz_p_negate(y.x, self.x)
        _sig_off
        return y

    def __cmp__(ntl_zz_p self, other):
        ## TODO: typecheck, think of what cmp should
        ## return if types are different. in the interim:
        if not PY_TYPE_CHECK(other, ntl_zz_p):
            return -1
        if (NTL_zz_p_DOUBLE_EQUALS(self.x, (<ntl_zz_p>other).x)):
            return 0
        else:
            return -1

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
        cdef ntl_zz_p y
        y = self._new()
        zz_p_sqr(y.x, self.x)
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
        ## TODO
        return zz_p_rep(self.x) == long(0)

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
        return zz_p_rep(self.x) == long(1)
        ## TODO
##        return 0
##        return zz_pX_IsOne(self.x)

    def clear(self):
        """
        Reset this element to 0.

        EXAMPLES:
        """
        ## TODO: check this
        zz_p_clear(self.x)


