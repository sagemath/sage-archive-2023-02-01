"""
Wrapper for NTL's polynomials over finite ring extensions of $\Z / p\Z.$

AUTHORS:
  -- David Roe (2007-10-10)
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                          David Roe     <roed@math.harvard.edu>
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
include 'misc.pxi'
include 'decl.pxi'

from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class
from sage.libs.ntl.ntl_ZZ_pEContext import ntl_ZZ_pEContext

from sage.libs.ntl.ntl_ZZ import unpickle_class_args

##############################################################################
#
# ZZ_pEX  -- polynomials over an extension of the integers modulo p
#
##############################################################################

cdef class ntl_ZZ_pEX:
    r"""
    The class \class{ZZ_pEX} implements polynomials over finite ring extensions of $\Z / p\Z$.

    It can be used, for example, for arithmentic in GF(p^n)[X].
    However, except where mathematically necessary (e.g., GCD computations),
    ZZ_pE need not be a field.
    """
    # See ntl_ZZ_pEX.pxd for definition of data members
    def __init__(self, v=None, modulus=None):
        """
        EXAMPLES:
            sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
            sage: a = ntl.ZZ_pE([3,2], c)
            sage: b = ntl.ZZ_pE([1,2], c)
            sage: f = ntl.ZZ_pEX([a, b, b])
            sage: f
            [[3 2] [1 2] [1 2]]
            sage: g = ZZ_pEX([0,0,0], c); g
            []
            sage: g[10]=5
            sage: g
            [[] [] [] [] [] [] [] [] [] [] [5]]
            sage: g[10]
            [5]
        """
        if modulus is None and v is None:
            raise ValueError, "You must specify a modulus when creating a ZZ_pEX."

        # self.c.restore_c()  ## Restoring the context is taken care of in __new__

        cdef ntl_ZZ_pE cc
        cdef Py_ssize_t i

        if v is None:
            return
        elif PY_TYPE_CHECK(v, list) or PY_TYPE_CHECK(v, tuple):
            for i from 0 <= i < len(v):
                x = v[i]
                if not PY_TYPE_CHECK(x, ntl_ZZ_pE):
                    cc = ntl_ZZ_pE(x,self.c)
                else:
                    if self.c is not (<ntl_ZZ_pE>x).c:
                        raise ValueError, "inconsistent moduli"
                    cc = x
                ZZ_pEX_SetCoeff(self.x, i, cc.x)
        else:
            raise NotImplementedError
            s = str(v).replace(',',' ').replace('L','')
            #_sig_on
            #ZZ_pEX_from_str(&self.x, s)
            #_sig_off

    def __new__(self, v=None, modulus=None):
        #################### WARNING ###################
        ## Before creating a ZZ_pEX, you must create a##
        ## ZZ_pEContext, and restore it.  In Python,  ##
        ## the error checking in __init__ will prevent##
        ## you from constructing an ntl_ZZ_pEX        ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = PY_NEW(ntl_ZZ_pEX) without    ##
        ## first restoring a ZZ_pEContext, which could##
        ## have unforetunate consequences.  See _new  ##
        ## defined below for an example of the right  ##
        ## way to short-circuit __init__ (or just call##
        ## _new in your own code).                    ##
        ################################################
        if modulus is None and v is None: # we also check for v is None so that a user can specify the modulus by v.
            ZZ_pEX_construct(&self.x)
            return
        if PY_TYPE_CHECK( modulus, ntl_ZZ_pEContext_class ):
            self.c = <ntl_ZZ_pEContext_class>modulus
        elif PY_TYPE_CHECK(v, ntl_ZZ_pEX):
            self.c = (<ntl_ZZ_pEX>v).c
        elif PY_TYPE_CHECK(v, ntl_ZZ_pE):
            self.c = (<ntl_ZZ_pE>v).c
        elif (PY_TYPE_CHECK(v, list) or PY_TYPE_CHECK(v, tuple)) and len(v) > 0:
            if PY_TYPE_CHECK(v[0], ntl_ZZ_pEX):
                self.c = (<ntl_ZZ_pEX>v[0]).c
            elif PY_TYPE_CHECK(v[0], ntl_ZZ_pE):
                self.c = (<ntl_ZZ_pEX>v[0]).c
            else:
                self.c = <ntl_ZZ_pEContext_class>ntl_ZZ_pEContext(modulus)
        else:
            self.c = <ntl_ZZ_pEContext_class>ntl_ZZ_pEContext(modulus)
        self.c.restore_c()
        ZZ_pEX_construct(&self.x)

    def __dealloc__(self):
        if <object>self.c is not None:
            self.c.restore_c()
        ZZ_pEX_destruct(&self.x)

    cdef ntl_ZZ_pEX _new(self):
        cdef ntl_ZZ_pEX r
        self.c.restore_c()
        r = PY_NEW(ntl_ZZ_pEX)
        r.c = self.c
        return r

    def __reduce__(self):
        """
        TEST:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: loads(dumps(f)) == f
        True
        """
        return make_ZZ_pEX, (self.list(), self.get_modulus_context())

    def __repr__(self):
        """
        Returns a string representation of self.

        TEST:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f
        [[3 2] [1 2] [1 2]]
        """
        self.c.restore_c()
        return ZZ_pEX_to_PyString(&self.x)
        #return string_delete(ZZ_pEX_to_str(&self.x))

    def __copy__(self):
        """
        Return a copy of self.

        TEST:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f
        [[3 2] [1 2] [1 2]]
        sage: y = f.copy()
        sage: y == f
        True
        sage: y is f
        False
        sage: f[0] = 0; y
        [[3 2] [1 2] [1 2]]
        """
        cdef ntl_ZZ_pEX r = self._new()
        #self.c.restore_c() ## _new() restores
        r.x = self.x
        return r

    def copy(self):
        """
        Return a copy of self.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f
        [[3 2] [1 2] [1 2]]
        sage: y = f.copy()
        sage: y == f
        True
        sage: y is f
        False
        sage: f[0] = 0; y
        [[3 2] [1 2] [1 2]]
        """
        return self.__copy__()

    def get_modulus_context(self):
        """
        Returns the structure that holds the underlying NTL modulus.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f.get_modulus_context()
        NTL modulus [1 1 1] (mod 7)
        """
        return self.c

    def __setitem__(self, long i, a):
        r"""
        Sets the ith coefficient of self to be a.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f[1] = 4; f
        [[3 2] [4] [1 2]]
        """
        if i < 0:
            raise IndexError, "index (i=%s) must be >= 0"%i
        cdef ntl_ZZ_pE _a
        if PY_TYPE_CHECK(a, ntl_ZZ_pE):
            _a = <ntl_ZZ_pE> a
        else:
            _a = ntl_ZZ_pE(a,self.c)
        self.c.restore_c()
        ZZ_pEX_SetCoeff(self.x, i, _a.x)

    def __getitem__(self, long i):
        r"""
        Returns the ith coefficient of self.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f[0]
        [3 2]
        """
        if i < 0:
            raise IndexError, "index (=%s) must be >= 0"%i
        cdef ntl_ZZ_pE r
        _sig_on
        self.c.restore_c()
        r = PY_NEW(ntl_ZZ_pE)
        r.c = self.c
        r.x = ZZ_pEX_coeff( self.x, i)
        _sig_off
        return r

    def list(self):
        """
        Return list of entries as a list of ntl_ZZ_pEs.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f.list()
        [[3 2], [1 2], [1 2]]
        """
        # This function could be sped up by using the list API and not restoring the context each time.
        # Or by using self.x.rep directly.
        self.c.restore_c()
        cdef Py_ssize_t i
        return [self[i] for i from 0 <= i <= self.degree()]

    def __add__(ntl_ZZ_pEX self, ntl_ZZ_pEX other):
        """
        Adds self and other.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: g = ntl.ZZ_pEX([-b, a])
        sage: f + g
        [[2] [4 4] [1 2]]
        """
        if self.c is not other.c:
            raise ValueError, "You can not perform arithmetic with elements of different moduli."
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        # self.c.restore_c() # _new restores the context
        ZZ_pEX_add(r.x, self.x, other.x)
        _sig_off
        return r

    def __sub__(ntl_ZZ_pEX self, ntl_ZZ_pEX other):
        """
        Subtracts other from self.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: g = ntl.ZZ_pEX([-b, a])
        sage: f - g
        [[4 4] [5] [1 2]]
        """
        if self.c is not other.c:
            raise ValueError, "You can not perform arithmetic with elements of different moduli."
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        # self.c.restore_c() # _new restores the context
        ZZ_pEX_sub(r.x, self.x, other.x)
        _sig_off
        return r

    def __mul__(ntl_ZZ_pEX self, ntl_ZZ_pEX other):
        """
        Returns the product self * other.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: g = ntl.ZZ_pEX([-b, a])
        sage: f * g
        [[1 3] [1 1] [2 4] [6 4]]
        sage: c2 = ntl.ZZ_pEContext(ntl.ZZ_pX([4,1,1], 5)) # we can mix up the moduli
        sage: x = c2.ZZ_pEX([2,4])
        sage: x^2
        [[4] [1] [1]]
        sage: f * g # back to the first one and the ntl modulus gets reset correctly
        [[1 3] [1 1] [2 4] [6 4]]
        """
        if self.c is not other.c:
            raise ValueError, "You can not perform arithmetic with elements of different moduli."
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        # self.c.restore_c() # _new() restores the context
        ZZ_pEX_mul(r.x, self.x, other.x)
        _sig_off
        return r

    def __div__(ntl_ZZ_pEX self, ntl_ZZ_pEX other):
        """
        Compute quotient self / other, if the quotient is a polynomial.
        Otherwise an Exception is raised.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a^2, -a*b-a*b, b^2])
        sage: g = ntl.ZZ_pEX([-a, b])
        sage: f / g
        [[4 5] [1 2]]
        sage: g / f
        Traceback (most recent call last):
        ...
        ArithmeticError: self (=[[4 5] [1 2]]) is not divisible by other (=[[5 1] [2 6] [4]])
        """
        if self.c is not other.c:
            raise ValueError, "You can not perform arithmetic with elements of different moduli."
        cdef int divisible
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        #self.c.restore_c() # _new restores context
        divisible = ZZ_pEX_divide(r.x, self.x, other.x)
        _sig_off
        if not divisible:
            raise ArithmeticError, "self (=%s) is not divisible by other (=%s)"%(self, other)
        return r

    def __mod__(ntl_ZZ_pEX self, ntl_ZZ_pEX other):
        """
        Given polynomials a, b in ZZ_pE[X], if p is prime and the defining modulus irreducible,
        there exist polynomials q, r in ZZ_pE[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns r.

        If p is not prime or the modulus is not irreducible, this function may raise a
        RuntimeError due to division by a noninvertible element of ZZ_p.

        EXAMPLES:
        sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([-5, 0, 1], 5^10))
        sage: a = c.ZZ_pE([5, 1])
        sage: b = c.ZZ_pE([4, 99])
        sage: f = c.ZZ_pEX([a, b])
        sage: g = c.ZZ_pEX([a^2, -b, a + b])
        sage: g % f
        [[1864280 2123186]]
        sage: f % g
        [[5 1] [4 99]]
        sage: f^2 % g
        ZZ_p: division by non-invertible element
        Traceback (most recent call last):
        ...
        RuntimeError:
        """
        if self.c is not other.c:
            raise ValueError, "You can not perform arithmetic with elements of different moduli."
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        # self.c.restore_c() # _new() restores the context
        ZZ_pEX_rem(r.x, self.x, other.x)
        _sig_off
        return r

    def quo_rem(self, ntl_ZZ_pEX other):
        """
        Given polynomials a, b in ZZ_pE[X], if p is prime and the defining modulus irreducible,
        there exist polynomials q, r in ZZ_pE[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns (q, r).

        If p is not prime or the modulus is not irreducible, this function may raise a
        RuntimeError due to division by a noninvertible element of ZZ_p.

        EXAMPLES:
        sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([-5, 0, 1], 5^10))
        sage: a = c.ZZ_pE([5, 1])
        sage: b = c.ZZ_pE([4, 99])
        sage: f = c.ZZ_pEX([a, b])
        sage: g = c.ZZ_pEX([a^2, -b, a + b])
        sage: g.quo_rem(f)
        ([[4947544 2492106] [4469276 6572944]], [[1864280 2123186]])
        sage: f.quo_rem(g)
        ([], [[5 1] [4 99]])
        sage: (f^2).quo_rem(g)
        ZZ_p: division by non-invertible element
        Traceback (most recent call last):
        ...
        RuntimeError:
        """
        if self.c is not other.c:
            raise ValueError, "You can not perform arithmetic with elements of different moduli."
        cdef ntl_ZZ_pEX r = self._new()
        cdef ntl_ZZ_pEX q = self._new()
        _sig_on
        # self.c.restore_c() # _new() restores the context
        ZZ_pEX_DivRem(q.x, r.x, self.x, other.x)
        _sig_off
        return q,r

    def square(self):
        """
        Return f^2.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f.square()
        [[5 1] [5 1] [2 1] [1] [4]]
        """
        # self.c.restore_c() # _new() restores the context
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_sqr(r.x, self.x)
        _sig_off
        return r

    def __pow__(ntl_ZZ_pEX self, long n, ignored):
        """
        Return the n-th nonnegative power of self.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f^5
        [[5 1] [2 6] [4 5] [5 1] [] [6 2] [2 3] [0 1] [1 4] [3 6] [2 4]]
        """
        self.c.restore_c()
        if n < 0:
            raise NotImplementedError
        import sage.rings.arith
        return sage.rings.arith.generic_power(self, n, ntl_ZZ_pEX([[1]],self.c))

    def __cmp__(ntl_ZZ_pEX self, ntl_ZZ_pEX other):
        """
        Decide whether or not self and other are equal.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: g = ntl.ZZ_pEX([a, b, b, 0])
        sage: f == g
        True
        sage: g = ntl.ZZ_pEX([a, b, a])
        sage: f == g
        False
        """
        self.c.restore_c()
        if ZZ_pEX_equal(self.x, other.x):
            return 0
        return -1

    def is_zero(self):
        """
        Return True exactly if this polynomial is 0.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f.is_zero()
        False
        sage: f = ntl.ZZ_pEX([0,0,7], c)
        sage: f.is_zero()
        True
        """
        self.c.restore_c()
        return bool(ZZ_pEX_IsZero(self.x))

    def is_one(self):
        """
        Return True exactly if this polynomial is 1.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f.is_one()
        False
        sage: f = ntl.ZZ_pEX([1, 0, 0], c)
        sage: f.is_one()
        True
        """
        self.c.restore_c()
        return bool(ZZ_pEX_IsOne(self.x))

    def is_monic(self):
        """
        Return True exactly if this polynomial is monic.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: f.is_monic()
        False
        sage: f = ntl.ZZ_pEX([a, b, 1], c)
        sage: f.is_monic()
        True
        """
        self.c.restore_c()
        # The following line is what we should have.  However, strangely this is *broken*
        # on PowerPC Intel in NTL, so we program around
        # the problem.  (William Stein)
        #return bool(ZZ_pEX_is_monic(self.x))

        if ZZ_pEX_IsZero(self.x):
             return False
        cdef ZZ_pE_c x = ZZ_pEX_LeadCoeff(self.x)
        return bool(ZZ_pE_IsOne(x))

    def __neg__(self):
        """
        Return the negative of self.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: -f
        [[4 5] [6 5] [6 5]]
        """
        cdef ntl_ZZ_pEX r = self._new()
        # self.c.restore_c() # _new() calls restore
        ZZ_pEX_negate(r.x, self.x)
        return r

    def left_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the left n positions.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b]); f
        [[3 2] [1 2] [1 2]]
        sage: f.left_shift(2)
        [[] [] [3 2] [1 2] [1 2]]
        sage: f.left_shift(5)
        [[] [] [] [] [] [3 2] [1 2] [1 2]]

        A negative left shift is a right shift.
        sage: f.left_shift(-2)
        [[1 2]]
        """
        # self.c.restore_c() # _new() calls restore
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_LeftShift(r.x, self.x, n)
        _sig_off
        return r

    def right_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the right n positions.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b]); f
        [[3 2] [1 2] [1 2]]
        sage: f.right_shift(2)
        [[1 2]]
        sage: f.right_shift(5)
        []

        A negative right shift is a left shift.
        sage: f.left_shift(-5)
        [[] [] [] [] [] [3 2] [1 2] [1 2]]
        """
        # self.c.restore_c() # _new() calls restore
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_RightShift(r.x, self.x, n)
        _sig_off
        return r

    def gcd(self, ntl_ZZ_pEX other, check=True):
        """
        Returns gcd(self, other) if we are working over a field.

        NOTE: Does not work if p is not prime or if the modulus is not irreducible.

        EXAMPLES:
        sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 11))
        sage: a = ntl.ZZ_pE([3,2], c)
        sage: b = ntl.ZZ_pE([1,2], c)
        sage: f = ntl.ZZ_pEX([a, b, b])
        sage: g = f^2
        sage: h = f^3
        sage: g.gcd(h)
        [[2 1] [8 1] [9 1] [2] [1]]
        sage: f^2
        [[5 8] [9 8] [6 8] [5] [8]]
        sage: eight = ntl.ZZ_pEX([[8]], c)
        sage: f^2 / eight
        [[2 1] [8 1] [9 1] [2] [1]]
        """
        #If check = True, need to check that ZZ_pE is a field.
        self.c.restore_c()
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_GCD(r.x, self.x, other.x)
        _sig_off
        return r

    def xgcd(self, ntl_ZZ_pEX other):
        """
        Returns r,s,t such that r = s*self + t*other.

        Here r is the resultant of a and b; if r != 0, then this
        function computes s and t such that: a*s + b*t = r; otherwise
        s and t are both 0.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = c.ZZ_pX([1,2,3]) * c.ZZ_pX([4,5])**2
            sage: g = c.ZZ_pX([1,1,1])**3 * c.ZZ_pX([1,2,3])
            sage: f.xgcd(g)   # nothing since they are not coprime
            ([6 12 1], [15 13 6 8 7 9], [4 13])

        In this example the input quadratic polynomials have a common root modulo 13.
            sage: f = c.ZZ_pX([5,0,1])
            sage: g = c.ZZ_pX([18,0,1])
            sage: f.xgcd(g)
            ([1], [13], [4])
            """
        self.c.restore_c()
        cdef ntl_ZZ_pEX s = self._new()
        cdef ntl_ZZ_pEX t = self._new()
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_XGCD(r.x, s.x, t.x, self.x, other.x)
        _sig_off
        return (r,s,t)

    def degree(self):
        """
        Return the degree of this polynomial.  The degree of the 0
        polynomial is -1.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([5,0,1])
            sage: f.degree()
            2
            sage: f = c.ZZ_pX(range(100))
            sage: f.degree()
            99
            sage: f = c.ZZ_pX()
            sage: f.degree()
            -1
            sage: f = c.ZZ_pX([1])
            sage: f.degree()
            0
        """
        self.c.restore_c()
        return ZZ_pEX_deg(self.x)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([3,6,9])
            sage: f.leading_coefficient()
            9
            sage: f = c.ZZ_pX()
            sage: f.leading_coefficient()
            0
        """
        self.c.restore_c()
        cdef long i
        i = ZZ_pEX_deg(self.x)
        return self[i]

    def constant_term(self):
        """
        Return the constant coefficient of this polynomial.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([3,6,9])
            sage: f.constant_term()
            3
            sage: f = c.ZZ_pX()
            sage: f.constant_term()
            0
        """
        self.c.restore_c()
        return self[0]

    def set_x(self):
        """
        Set this polynomial to the monomial "x".

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX()
            sage: f.set_x()
            sage: f
            [0 1]
            sage: g = c.ZZ_pX([0,1])
            sage: f == g
            True

        Though f and g are equal, they are not the same objects in memory:
            sage: f is g
            False
        """
        self.c.restore_c()
        ZZ_pEX_SetX(self.x)

    def is_x(self):
        """
        True if this is the polynomial "x".

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX()
            sage: f.set_x()
            sage: f.is_x()
            True
            sage: f = c.ZZ_pX([0,1])
            sage: f.is_x()
            True
            sage: f = c.ZZ_pX([1])
            sage: f.is_x()
            False
        """
        return bool(ZZ_pEX_IsX(self.x))

    def derivative(self):
        """
        Return the derivative of this polynomial.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([1,7,0,13])
            sage: f.derivative()
            [7 0 19]
        """
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_diff(r.x, self.x)
        _sig_off
        return r

    #def factor(self, verbose=False):
    #    cdef ZZ_pX_c** v
    #    cdef long* e
    #    cdef long i, n
    #    _sig_on
    #    ZZ_pX_factor(&v, &e, &n, &self.x, verbose)
    #    _sig_off
    #    F = []
    #    for i from 0 <= i < n:
    #        F.append((make_ZZ_pX(v[i], self.c), e[i]))
    #    free(v)
    #    free(e)
    #    return F

    #def linear_roots(self):
    #    """
    #    Assumes that input is monic, and has deg(f) distinct roots.
    #    Returns the list of roots.
    #    """
    #    cdef ZZ_p_c** v
    #    cdef long i, n
    #    _sig_on
    #    ZZ_pX_linear_roots(&v, &n, &self.x)
    #    _sig_off
    #    F = []
    #    for i from 0 <= i < n:
    #        F.append(make_ZZ_p(v[i], self.c))
    #    free(v)
    #    return F

    def reverse(self, hi=None):
        """
        Return the polynomial obtained by reversing the coefficients
        of this polynomial.  If hi is set then this function behaves
        as if this polynomial has degree hi.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([1,2,3,4,5])
            sage: f.reverse()
            [5 4 3 2 1]
            sage: f.reverse(hi=10)
            [0 0 0 0 0 0 5 4 3 2 1]
            sage: f.reverse(hi=2)
            [3 2 1]
            sage: f.reverse(hi=-2)
            []
        """
        cdef ntl_ZZ_pEX r = self._new()
        if not (hi is None):
            ZZ_pEX_reverse_hi(r.x, self.x, int(hi))
        else:
            ZZ_pEX_reverse(r.x, self.x)
        return r

    def truncate(self, long m):
        """
        Return the truncation of this polynomial obtained by
        removing all terms of degree >= m.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([1,2,3,4,5])
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
        cdef ntl_ZZ_pEX r = self._new()
        if m > 0:
            _sig_on
            ZZ_pEX_trunc(r.x, self.x, m)
            _sig_off
        return r

    def multiply_and_truncate(self, ntl_ZZ_pEX other, long m):
        """
        Return self*other but with terms of degree >= m removed.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([1,2,3,4,5])
            sage: g = c.ZZ_pX([10])
            sage: f.multiply_and_truncate(g, 2)
            [10]
            sage: g.multiply_and_truncate(f, 2)
            [10]
        """
        cdef ntl_ZZ_pEX r = self._new()
        if m > 0:
            _sig_on
            ZZ_pEX_MulTrunc(r.x, self.x, other.x, m)
            _sig_off
        return r

    def square_and_truncate(self, long m):
        """
        Return self*self but with terms of degree >= m removed.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([1,2,3,4,5])
            sage: f.square_and_truncate(4)
            [1 4 10]
            sage: (f*f).truncate(4)
            [1 4 10]
        """
        cdef ntl_ZZ_pEX r = self._new()
        if m > 0:
            _sig_on
            ZZ_pEX_SqrTrunc(r.x, self.x, m)
            _sig_off
        return r

    def invert_and_truncate(self, long m):
        """
        Compute and return the inverse of self modulo $x^m$.
        The constant term of self must be invertible.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([1,2,3,4,5,6,7])
            sage: f.invert_and_truncate(20)
            [1 18 1 0 0 0 0 8 17 2 13 0 0 0 4 0 17 10 9]
            sage: g = f.invert_and_truncate(20)
            sage: g * f
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 4 4 3]
        """
        if m < 0:
            raise ArithmeticError, "m (=%s) must be positive"%m
        #Need to check here if constant term is invertible
        cdef ntl_ZZ_pEX r = self._new()
        if m > 0:
            _sig_on
            ZZ_pEX_InvTrunc(r.x, self.x, m)
            _sig_off
        return r

    def multiply_mod(self, ntl_ZZ_pEX other, ntl_ZZ_pEX modulus):
        """
        Return self*other % modulus.  The modulus must be monic with
        deg(modulus) > 0, and self and other must have smaller degree.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: modulus = c.ZZ_pX([1,2,0,1])    # must be monic
            sage: g = c.ZZ_pX([-1,0,1])
            sage: h = c.ZZ_pX([3,7,13])
            sage: h.multiply_mod(g, modulus)
            [10 6 4]
        """
        self.c.restore_c()
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_MulMod(r.x, self.x, other.x, modulus.x)
        _sig_off
        return r

    def trace_mod(self, ntl_ZZ_pEX modulus):
        """
        Return the trace of this polynomial modulo the modulus.
        The modulus must be monic, and of positive degree degree bigger
        than the the degree of self.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = c.ZZ_pX([1,2,0,3])
            sage: mod = c.ZZ_pX([5,3,-1,1,1])
            sage: f.trace_mod(mod)
            3
        """
        self.c.restore_c()
        cdef ntl_ZZ_pE r = ntl_ZZ_pE(modulus = self.c)
        _sig_on
        ZZ_pEX_TraceMod(r.x, self.x, modulus.x)
        _sig_off
        return r

    #def trace_list(self):
    #    """
    #    Return the list of traces of the powers $x^i$ of the
    #    monomial x modulo this polynomial for i = 0, ..., deg(f)-1.
    #    This polynomial must be monic.
    #
    #    EXAMPLES:
    #        sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
    #        sage: f = c.ZZ_pX([1,2,0,3,0,1])
    #        sage: f.trace_list()
    #        [5, 0, 14, 0, 10]
    #
    #        The input polynomial must be monic or a ValueError is raised:
    #        sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
    #        sage: f = c.ZZ_pX([1,2,0,3,0,2]
    #        sage: f.trace_list()
    #        Traceback (most recent call last):
    #        ...
    #        ValueError: polynomial must be monic.
    #    """
    #    self.c.restore_c()
    #    if not self.is_monic():
    #        raise ValueError, "polynomial must be monic."
    #    cdef long N = self.degree()
    #    cdef vec_ZZ_pE_c
    #    _sig_on
    #    cdef char* t
    #    t = ZZ_pX_trace_list(&self.x)
    #    return eval(string(t).replace(' ', ','))

    def resultant(self, ntl_ZZ_pEX other):
        """
        Return the resultant of self and other.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = c.ZZ_pX([17,0,1,1])
            sage: g = c.ZZ_pX([34,-17,18,2])
            sage: f.resultant(g)
            0
        """
        self.c.restore_c()
        cdef ntl_ZZ_pE r = ntl_ZZ_pE(modulus = self.c)
        _sig_on
        ZZ_pEX_resultant(r.x, self.x, other.x)
        _sig_off
        return r

    def norm_mod(self, ntl_ZZ_pEX modulus):
        """
        Return the norm of this polynomial modulo the modulus.  The
        modulus must be monic, and of positive degree strictly greater
        than the degree of self.

        EXAMPLE:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = c.ZZ_pX([1,2,0,3])
            sage: mod = c.ZZ_pX([-5,2,0,0,1])
            sage: f.norm_mod(mod)
            11

        The norm is the constant term of the characteristic polynomial.
            sage: f.charpoly_mod(mod)
            [11 1 8 14 1]
        """
        self.c.restore_c()
        cdef ntl_ZZ_pE r = ntl_ZZ_pE(modulus = self.c)
        _sig_on
        ZZ_pEX_NormMod(r.x, self.x, modulus.x)
        _sig_off
        return r

    def discriminant(self):
        r"""
        Return the discriminant of a=self, which is by definition
        $$
                (-1)^{m(m-1)/2} {\mbox{\tt resultant}}(a, a')/lc(a),
        $$
        where m = deg(a), and lc(a) is the leading coefficient of a.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = c.ZZ_pX([1,2,0,3])
            sage: f.discriminant()
            1
        """
        self.c.restore_c()
        cdef long m

        c = ~self.leading_coefficient()
        m = self.degree()
        if (m*(m-1)/2) % 2:
            c = -c
        return c*self.resultant(self.derivative())

    #def charpoly_mod(self, ntl_ZZ_pEX modulus):
    #    """
    #    Return the characteristic polynomial of this polynomial modulo
    #    the modulus.  The modulus must be monic of degree bigger than
    #    self.
    #
    #    EXAMPLES:
    #        sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
    #        sage: f = c.ZZ_pX([1,2,0,3])
    #        sage: mod = c.ZZ_pX([-5,2,0,0,1])
    #        sage: f.charpoly_mod(mod)
    #        [11 1 8 14 1]
    #    """
    #    self.c.restore_c()
    #    cdef ntl_ZZ_pEX r = self._new()
    #    _sig_on
    #    ZZ_pEX_charpoly_mod(&r.x, &self.x, &modulus.x)
    #    _sig_off
    #    return r

    def minpoly_mod(self, ntl_ZZ_pEX modulus):
        """
        Return the minimal polynomial of this polynomial modulo the
        modulus.  The modulus must be monic of degree bigger than
        self.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = c.ZZ_pX([0,0,1])
            sage: g = f*f
            sage: f.charpoly_mod(g)
            [0 0 0 0 1]

        However, since $f^2 = 0$ moduluo $g$, its minimal polynomial
        is of degree $2$.
            sage: f.minpoly_mod(g)
            [0 0 1]
        """
        self.c.restore_c()
        cdef ntl_ZZ_pEX r = self._new()
        _sig_on
        ZZ_pEX_MinPolyMod(r.x, self.x, modulus.x)
        _sig_off
        return r

    def clear(self):
        """
        Reset this polynomial to 0.  Changes this polynomial in place.

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = c.ZZ_pX([1,2,3])
            sage: f
            [1 2 3]
            sage: f.clear()
            sage: f
            []
        """
        self.c.restore_c()
        _sig_on
        ZZ_pEX_clear(self.x)
        _sig_off

    def preallocate_space(self, long n):
        """
        Pre-allocate spaces for n coefficients.  The polynomial that f
        represents is unchanged.  This is useful if you know you'll be
        setting coefficients up to n, so memory isn't re-allocated as
        the polynomial grows.  (You might save a millisecond with this
        function.)

        EXAMPLES:
            sage: c=ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = c.ZZ_pX([1,2,3])
            sage: f.preallocate_space(20)
            sage: f
            [1 2 3]
            sage: f[10]=5  # no new memory is allocated
            sage: f
            [1 2 3 0 0 0 0 0 0 0 5]
        """
        self.c.restore_c()
        _sig_on
        self.x.SetMaxLength(n)
        _sig_off


def make_ZZ_pEX(v, modulus):
    """
    Here for unpickling.

    EXAMPLES:
    sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([-5,0,1],25))
    sage: sage.libs.ntl.ntl_ZZ_pE.make_ZZ_pE(ntl.ZZ_pE, [4,3], c)
    [4 3]
    sage: type(sage.libs.ntl.ntl_ZZ_pE.make_ZZ_pE(ntl.ZZ_pE, [4,3], c))
    <type 'sage.libs.ntl.ntl_ZZ_pE.ntl_ZZ_pE'>
    """

    return ntl_ZZ_pEX(v, modulus)
