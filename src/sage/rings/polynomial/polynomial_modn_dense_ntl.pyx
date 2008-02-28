"""
Dense univariate polynomials over Z/nZ, implemented using NTL.

AUTHORS:
    -- Robert Bradshaw: Split off from polynomial_element_generic.py (2007-09)
    -- Robert Bradshaw: Major rewrite to use NTL directly (2007-09)

"""

################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

from sage.rings.polynomial.polynomial_element import is_Polynomial, Polynomial_generic_dense

from sage.libs.all import pari, pari_gen

from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZX, zero_ZZX, ZZ_p, ZZ_pX ##, set_modulus
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.integer_mod import IntegerMod_abstract

from sage.rings.fraction_field_element import FractionFieldElement
import sage.rings.polynomial.polynomial_ring

from sage.rings.infinity import infinity

import polynomial_singular_interface
from sage.interfaces.all import singular as singular_default

from sage.structure.element import generic_power, canonical_coercion, bin_op

def make_element(parent, args):
    return parent(*args)

include "../../ext/stdsage.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/cdefs.pxi"

zz_p_max = NTL_SP_BOUND

cdef class Polynomial_dense_mod_n(Polynomial):
    """
    A dense polynomial over the integers modulo n, where n is composite.

    Much of the underlying arithmetic is done using NTL.

    EXAMPLES:

        sage: R.<x> = PolynomialRing(Integers(16))
        sage: f = x^3 - x + 17
        sage: f^2
        x^6 + 14*x^4 + 2*x^3 + x^2 + 14*x + 1

        sage: loads(f.dumps()) == f
        True

        sage: R.<x> = Integers(100)[]
        sage: p = 3*x
        sage: q = 7*x
        sage: p+q
        10*x
        sage: R.<x> = Integers(8)[]
        sage: parent(p)
        Univariate Polynomial Ring in x over Ring of integers modulo 100
        sage: p + q
        10*x

    """
    def __init__(self, parent, x=None, check=True,
                 is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        cdef Polynomial_dense_mod_n numer, denom

        if construct:
            if isinstance(x, ZZ_pX):
                self.__poly = x
                return
            self.__poly = ZZ_pX(x, parent.modulus())
            return

        self.__poly = ZZ_pX([], parent.modulus())

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = (<Polynomial_dense_modn_ntl_zz>x).__poly.__copy__()
                return
            else:
                R = parent.base_ring()
                x = [ZZ(R(a)) for a in x.list()]
                check = False

        elif isinstance(x, dict):
            x = self._dict_to_list(x, R(0))


        elif isinstance(x, ZZX):
            self.__poly = x.copy()
            return

        elif isinstance(x, pari_gen):
            x = [ZZ(w) for w in x.Vecrev()]
            check = False

        elif isinstance(x, FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_dense_mod_n):
            if x.denominator().is_unit():
                numer = x.numerator()
                denom = x.denominator().inverse_of_unit()
                x = numer.__poly * denom.__poly
                check = False
            else:
                raise TypeError, "Denominator not a unit."

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            R = parent.base_ring()
            x = [ZZ(R(a)) for a in x]

        self.__poly = ZZ_pX(x, parent.modulus())

##    def _ntl_set_modulus(self):
##        self.parent()._ntl_set_modulus()

    def __reduce__(self):
        return make_element, (self.parent(), (self.list(), False, self.is_gen()))

    def int_list(self):
        return eval(str(self.__poly).replace(' ',','))

    def _pari_(self, variable=None):
        """
        EXAMPLES:
            sage: t = PolynomialRing(IntegerModRing(17),"t").gen()
            sage: f = t^3 + 3*t - 17
            sage: pari(f)
            Mod(1, 17)*t^3 + Mod(3, 17)*t
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.int_list()).Polrev(variable) * \
               pari(1).Mod(self.parent().base_ring().order())

    def ntl_ZZ_pX(self):
        r"""
        Return underlying NTL representation of this polynomial.
        Additional ``bonus'' functionality is available through this
        function.

        WARNING:
        You must call \code{ntl.set_modulus(ntl.ZZ(n))} before doing
        arithmetic with this object!
        """
        return self.__poly

    def __getitem__(self, n):
        return self.parent().base_ring()(self.__poly[n]._sage_())

    def __getslice__(self, i, j):
        R = self.base_ring()
        if i < 0:
            i = 0
        if j > self.__poly.degree()+1:
            j = self.__poly.degree()+1
        v = [R(self.__poly[k]._sage_()) for k in range(i,j)]
        return self.parent()([0]*int(i) + v)

    def _unsafe_mutate(self, n, value):
        n = int(n)
        if n < 0:
            raise IndexError, "n must be >= 0"
##        self._ntl_set_modulus()
        self.__poly[n] = int(value)

    def _pow(self, n):
        n = int(n)
##        self._ntl_set_modulus()
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
##        self._ntl_set_modulus()
        return self.parent()(self.__poly + (<Polynomial_dense_mod_n>right).__poly, construct=True)

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: x = PolynomialRing(Integers(100), 'x').0
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 + 90*x^2 + 32*x + 68
        """
##        self._ntl_set_modulus()
        return self.parent()(self.__poly * (<Polynomial_dense_mod_n>right).__poly, construct=True)

    def _rmul_(self, c):
        try:
##            self._ntl_set_modulus()
            return self.parent()(ZZ_pX([c], self.parent().modulus()) * self.__poly, construct=True)
        except RuntimeError, msg: # should this realy be a TypeError
            raise TypeError, msg

    def _lmul_(self, c):
        try:
##            self._ntl_set_modulus()
            return self.parent()(ZZ_pX([c], self.parent().modulus()) * self.__poly, construct=True)
        except RuntimeError, msg: # should this realy be a TypeError
            raise TypeError, msg

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_dense_mod_n):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
##        self._ntl_set_modulus()
        v = self.__poly.quo_rem((<Polynomial_dense_mod_n>right).__poly)
        P = self.parent()
        return (P(v[0], construct=True), P(v[1], construct=True) )

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power $x^n$. If $n$ is negative,
        terms below $x^n$ will be discarded. Does not change this polynomial.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(Integers(12345678901234567890))
            sage: p = x^2 + 2*x + 4
            sage: p.shift(0)
             x^2 + 2*x + 4
            sage: p.shift(-1)
             x + 2
            sage: p.shift(-5)
             0
            sage: p.shift(2)
             x^4 + 2*x^3 + 4*x^2

        AUTHOR:
            -- David Harvey (2006-08-06)
        """
        if n == 0:
            return self
##        self._ntl_set_modulus()
        return self.parent()(self.__poly.left_shift(n),
                             construct=True)

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
##        self._ntl_set_modulus()
        return self.parent()(self.__poly - (<Polynomial_dense_mod_n>right).__poly, construct=True)

    def __floordiv__(self, right):
        if is_Polynomial(right) and right.is_constant() and \
                         right[0] in self.parent().base_ring():
            d = right[0]
        elif (right in self.parent().base_ring()):
            d = right
        else:
            return Polynomial.__floordiv__(self, right)
        return self.parent()([c // d for c in self.list()], construct=True)

##     def __copy__(self):
##         self.parent()._ntl_set_modulus()
##         f = self.parent()()
##         f.__poly = self.__poly.copy()
##         return f

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.degree(), -1)

    def is_irreducible(self):
        return bool(self._pari_().polisirreducible())

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES:
            sage: _.<x> = Integers(100)[]
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [83, 3, 0, 1]
        """
        R = self.base_ring()
        return [R(x) for x in self.int_list()]

    def ntl_set_directly(self, v):
        r"""
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers modulo n are stored internally
        using NTL's ZZ_pX class.  Use this function to set the value
        of this polynomial using the NTL constructor, which is
        potentially \emph{very} fast.  The input v is either a vector
        of ints or a string of the form '[ n1 n2 n3 ... ]' where the
        ni are integers and there are no commas between them.  The
        optimal input format is the string format, since that's what
        NTL uses by default.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(Integers(100))
            sage: R([1,-2,3])
            3*x^2 + 98*x + 1
            sage: f = R(0)
            sage: f.ntl_set_directly([1,-2,3])
            sage: f
            3*x^2 + 98*x + 1
            sage: f.ntl_set_directly('[1 -2 3 4]')
            sage: f
            4*x^3 + 3*x^2 + 98*x + 1
        """
        if self.is_gen():
            raise TypeError, "Cannot change the value of the generator."
##        self._ntl_set_modulus()
        self.__poly = ZZ_pX(v, self.parent().modulus())

    # Polynomial_singular_repr stuff, copied due to lack of multiple inheritance

    def _singular_(self, singular=singular_default, have_ring=False, force=False):
        if not have_ring:
            self.parent()._singular_(singular,force=force).set_ring() #this is expensive
        if self.__singular is not None:
            try:
                self.__singular._check_valid()
                if self.__singular.parent() is singular:
                    return self.__singular
            except (AttributeError, ValueError):
                pass
        return self._singular_init_(singular, have_ring=have_ring)

    def _singular_init_(self, singular=singular_default, have_ring=False, force=False):
        if not have_ring:
            self.parent()._singular_(singular,force=force).set_ring() #this is expensive
        self.__singular = singular(str(self))
        return self.__singular

    def lcm(self, singular=singular_default, have_ring=False):
        return polynomial_singular_interface.lcm_func(self, singular, have_ring)

    def resultant(self, other, variable=None):
        return polynomial_singular_interface.resultant_func(self, other, variable)




cdef class Polynomial_dense_modn_ntl_zz(Polynomial_dense_mod_n):

    def __init__(self, parent, v=None, check=True, is_gen=False, construct=False):
        r"""
        EXAMPLES:
            sage: R = Integers(5**21) ; S.<x> = R[]
            sage: S(1/4)
            357627868652344
        """
        if isinstance(v, Polynomial):
            if (<Element>v)._parent == parent:
                Polynomial.__init__(self, parent, is_gen=is_gen)
                self.x = (<Polynomial_dense_modn_ntl_zz>v).x
                self.c = (<Polynomial_dense_modn_ntl_zz>v).c
                return

        Polynomial_dense_mod_n.__init__(self, parent, v, check=check, is_gen=is_gen, construct=construct)
#        if check:
#            R = parent.base_ring()
#            v = [a if isinstance(a, (int, long, Integer, IntegerMod_abstract)) else R(a) for a in v]
        v = [a for a in self.__poly.list()]
        self.__poly = None # this will eventually go away
        cdef ntl_zz_pX ntl = ntl_zz_pX(v, parent.modulus()) # let it handle the hard work
        self.x = ntl.x
        self.c = ntl.c

    def __dealloc__(self):
        if <object>self.c is not None:
            self.c.restore_c()
        zz_pX_destruct(&self.x)

    def ntl_set_directly(self, v):
        # TODO: Get rid of this
        Polynomial_dense_mod_n.ntl_set_directly(self, v)
        # verbatim from __init__
        v = [int(a) for a in self.__poly.list()]
        self.__poly = None # this will eventually go away
        cdef ntl_zz_pX ntl = ntl_zz_pX(v, self._parent.modulus()) # let it handle the hard work
        self.x = ntl.x
        self.c = ntl.c

    cdef Polynomial_dense_modn_ntl_zz _new(self):
        cdef Polynomial_dense_modn_ntl_zz y = <Polynomial_dense_modn_ntl_zz>PY_NEW(Polynomial_dense_modn_ntl_zz)
        y.c = self.c
        y._parent = self._parent
        return y

    def int_list(self):
        """
        Returns the coefficents of self as efficiently as possible as a
        list of python ints.

        EXAMPLES:
            sage: R.<x> = Integers(100)[]
            sage: f = x^3 + 5
            sage: f.int_list()
            [5, 0, 0, 1]
            sage: [type(a) for a in f.int_list()]
            [<type 'int'>, <type 'int'>, <type 'int'>, <type 'int'>]
        """
        cdef long i
        return [ zz_p_rep(zz_pX_GetCoeff(self.x, i)) for i from 0 <= i <= zz_pX_deg(self.x) ]

    def __getitem__(self, n):
        """
        EXAMPLES:
            sage: R.<x> = Integers(100)[]
            sage: f = (x+2)^7
            sage: f[3]
            60
        """
        R = self._parent._base
        if n < 0 or n > zz_pX_deg(self.x):
            return R(0)
        else:
            return R(zz_p_rep(zz_pX_GetCoeff(self.x, n)))

    def _unsafe_mutate(self, n, value):
        self.c.restore_c()
        zz_pX_SetCoeff_long(self.x, n, value)

    def __getslice__(self, i, j):
        """
        Returns the monomials of self of degree i <= n < j.

        EXAMPLES:
            sage: R.<x> = Integers(100)[]
            sage: f = (x+2)^7
            sage: f[3:6]
            84*x^5 + 80*x^4 + 60*x^3
            sage: f[-5:50] == f
            True
            sage: f[6:]
            x^7 + 14*x^6
        """
        R = self.base_ring()
        if i < 0:
            i = 0
        if j > zz_pX_deg(self.x)+1:
            j = zz_pX_deg(self.x)+1
        v = [ zz_p_rep(zz_pX_GetCoeff(self.x, t)) for t from i <= t < j ]
        return Polynomial_dense_modn_ntl_zz(self._parent, v, check=False) << i

    cdef ModuleElement _add_c_impl(self, ModuleElement _right):
        """
        TESTS:
            sage: R.<x> = Integers(100)[]
            sage: (x+5) + (x^2 - 6)
            x^2 + x + 99
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 1000000
        if do_sig: _sig_on
        self.c.restore_c()
        zz_pX_add(r.x, self.x, right.x)
        if do_sig: _sig_off
        return r

    cdef ModuleElement _sub_c_impl(self, ModuleElement _right):
        """
        TESTS:
            sage: R.<x> = Integers(100)[]
            sage: (x+5) - (x^2 - 6)
            99*x^2 + x + 11
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 1000000
        if do_sig: _sig_on
        self.c.restore_c()
        zz_pX_sub(r.x, self.x, right.x)
        if do_sig: _sig_off
        return r

    cdef RingElement _mul_c_impl(self, RingElement _right):
        """
        TESTS:
            sage: R.<x> = Integers(100)[]
            sage: (x+5) * (x^2 - 1)
            x^3 + 5*x^2 + 99*x + 95
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 10000
        if do_sig: _sig_on
        self.c.restore_c()
        if self is right:
            zz_pX_sqr(r.x, self.x)
        else:
            zz_pX_mul(r.x, self.x, right.x)
        if do_sig: _sig_off
        return r

    cdef ModuleElement _rmul_c_impl(self, RingElement c):
        """
        TESTS:
            sage: R.<x> = Integers(100)[]
            sage: (x+5) * 3
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) > 100000
        if do_sig: _sig_on
        self.c.restore_c()
        zz_pX_rmul(r.x, self.x, c)
        if do_sig: _sig_off
        return r

    cdef ModuleElement _lmul_c_impl(self, RingElement c):
        """
        TESTS:
            sage: R.<x> = Integers(100)[]
            sage: 3 * (x+5)
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) > 100000
        if do_sig: _sig_on
        self.c.restore_c()
        zz_pX_lmul(r.x, c, self.x)
        if do_sig: _sig_off
        return r

    def __pow__(Polynomial_dense_modn_ntl_zz self, ee, dummy):
        """
        TESTS:
            sage: R.<x> = Integers(100)[]
            sage: (x-1)^5
            x^5 + 95*x^4 + 10*x^3 + 90*x^2 + 5*x + 99

            sage: R.<x> = Integers(101)[]
            sage: (x-1)^(-5)
            1/(x^5 + 96*x^4 + 10*x^3 + 91*x^2 + 5*x + 100)
        """
        cdef bint recip = 0, do_sig
        cdef long e = ee
        if e != ee:
            raise TypeError, "Only integral powers defined."
        elif e < 0:
            recip = 1
            e = -e
        if not self:
            if e == 0:
                raise ArithmeticError, "0^0 is undefined."
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        self.c.restore_c()
        if zz_pX_IsX(self.x):
            zz_pX_LeftShift(r.x, self.x, e-1)
        else:
            do_sig = zz_pX_deg(self.x) *e > 1000
            if do_sig: _sig_on
            zz_pX_power(r.x, self.x, e)
            if do_sig: _sig_off
        if recip:
            return ~r
        else:
            return r

    def quo_rem(self, right):
        """
        Returns $q$ and $r$, with the degree of $r$ less than the degree of $right$,
        such that $q * right + r = self$.

        EXAMPLES:
            sage: R.<x> = Integers(125)[]
            sage: f = x^5+1; g = (x+1)^2
            sage: q, r = f.quo_rem(g)
            sage: q
            x^3 + 123*x^2 + 3*x + 121
            sage: r
            5*x + 5
            sage: q*g + r
            x^5 + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or self._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self.quo_rem(right)
        cdef Polynomial_dense_modn_ntl_zz q = self._new()
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(denom.x) > 1000
        if do_sig: _sig_on
        self.c.restore_c()
        zz_pX_divrem(q.x, r.x, self.x, denom.x)
        if do_sig: _sig_off
        return q, r

    def __floordiv__(self, right):
        """
        Returns the whole part of self/right, without remainder.

        For q = n // d, we have deg(n - q*d) < deg(d)

        EXAMPLES:
            sage: R.<x> = Integers(25)[]
            sage: f = x^7 + 1; g = x^2 - 1
            sage: q = f // g; q
            x^5 + x^3 + x
            sage: f - q*g
            x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self // right
        cdef Polynomial_dense_modn_ntl_zz numer = <Polynomial_dense_modn_ntl_zz>self
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        cdef Polynomial_dense_modn_ntl_zz q = numer._new()
        cdef bint do_sig = zz_pX_deg(numer.x) + zz_pX_deg(denom.x) > 1000
        if do_sig: _sig_on
        numer.c.restore_c()
        zz_pX_div(q.x, numer.x, denom.x)
        if do_sig: _sig_off
        return q

    def __mod__(self, right):
        """
        EXAMPLES:
            sage: R.<x> = Integers(81)[]
            sage: f = x^7 + x + 1; g = x^3
            sage: r = f % g; r
            x + 1
            sage: g * x^4 + r
            x^7 + x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self % right
        cdef Polynomial_dense_modn_ntl_zz numer = <Polynomial_dense_modn_ntl_zz>self
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        cdef Polynomial_dense_modn_ntl_zz r = numer._new()
        cdef bint do_sig = zz_pX_deg(numer.x) + zz_pX_deg(denom.x) > 1000
        if do_sig: _sig_on
        numer.c.restore_c()
        zz_pX_mod(r.x, numer.x, denom.x)
        if do_sig: _sig_off
        return r

    def shift(self, n):
        """
        Shift self to left by $n$, which is multiplication by $x^n$,
        truncating if $n$ is negative.

        EXAMPLES:
            sage: R.<x> = Integers(77)[]
            sage: f = x^7 + x + 1
            sage: f.shift(1)
            x^8 + x^2 + x
            sage: f.shift(-1)
            x^6 + 1
            sage: f.shift(10).shift(-10) == f
            True
        """
        return self << n

    def __lshift__(Polynomial_dense_modn_ntl_zz self, long n):
        """
        TEST:
            sage: R.<x> = Integers(77)[]
            sage: f = x^5 + 2*x + 1
            sage: f << 3
            x^8 + 2*x^4 + x^3
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_LeftShift(r.x, self.x, n)
        return r

    def __rshift__(Polynomial_dense_modn_ntl_zz self, long n):
        """
        TEST:
            sage: R.<x> = Integers(77)[]
            sage: f = x^5 + 2*x + 1
            sage: f >> 3
            x^2
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_RightShift(r.x, self.x, n)
        return r

    def _derivative(self, var=None):
        r"""
        Returns the formal derivative of self with respect to var.

        var must be either the generator of the polynomial ring to which
        this polynomial belongs, or None (either way the behaviour is the
        same).

        SEE ALSO:
            self.derivative()

        EXAMPLES:
            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f._derivative()
            4*x^3 + 76
            sage: f._derivative(None)
            4*x^3 + 76

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

            sage: y = var("y")
            sage: f._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var is not self._parent.gen():
            raise ValueError, "cannot differentiate with respect to %s" % var

        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_diff(r.x, self.x)
        return r


    def reverse(self):
        """
        Reverses the coeffients of self. The reverse of f(x) is x^n f(1/x).

        The degree will go down if the constant term is zero.

        EXAMPLES:
            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f.reverse()
            76*x^4 + 76*x^3 + 1
            sage: f = x^3 - x
            sage: f.reverse()
            76*x^2 + 1
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_reverse(r.x, self.x)
        return r

    def is_gen(self):
        return zz_pX_IsX(self.x)

    def __nonzero__(self):
        """
        TESTS:
            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: not f
            False
            sage: not (x-x)
            True
        """
        return not zz_pX_IsZero(self.x)

    def valuation(self):
        """
        Returns the valuation of self, that is, the power of the
        lowest non-zero monomial of self.

        EXAMPLES:
            sage: R.<x> = Integers(10)[]
            sage: x.valuation()
            1
            sage: f = x-3; f.valuation()
            0
            sage: f = x^99; f.valuation()
            99
            sage: f = x-x; f.valuation()
            +Infinity
        """
        cdef long n
        for n from 0 <= n <= zz_pX_deg(self.x):
            if zz_p_rep(zz_pX_GetCoeff(self.x, n)):
                return n
        return infinity

    def degree(self):
        """
        EXAMPLES:
            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f.degree()
            4
            sage: f = 77*x + 1
            sage: f.degree()
            0
        """
        return zz_pX_deg(self.x)

    def truncate(self, long n):
        """
        Returns this polynomial mod $x^n$.

        EXAMPLES:
            sage: R.<x> = Integers(77)[]
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_trunc(r.x, self.x, n)
        return r

    def __call__(self, *args, **kwds):
        """
        Evaluate self at x. If x is a single argument coercable into
        the basering of self, this is done directly in NTL, otherwise
        the generic Polynomial call code is used.

        EXAMPLES:
            sage: R.<x> = Integers(100)[]
            sage: f = x^3+7
            sage: f(5)
            32
            sage: f(5r)
            32
            sage: f(mod(5, 1000))
            32
            sage: f(x)
            x^3 + 7
            sage: S.<y> = Integers(5)[]
            sage: f(y)
            y^3 + 2
        """
        if len(args) != 1 or len(kwds) != 0:
            return Polynomial.__call__(self, *args, **kwds)
        arg = args[0]
#        cdef zz_p_c x
        cdef ntl_zz_p fx = ntl_zz_p(0, self.c), x = None
        if PY_TYPE_CHECK(arg, int):
            x = ntl_zz_p(arg, self.c)
        elif PY_TYPE_CHECK(arg, Integer):
            x = ntl_zz_p(arg, self.c)
        elif PY_TYPE_CHECK(arg, Element):
            if <void *>self._parent._base == <void *>(<Element>arg)._parent: # c++ pointer hack
                x = ntl_zz_p(arg, self.c)
            else:
                map = self._parent._base.coerce_map_from((<Element>arg)._parent)
                if map is not None:
                    x = ntl_zz_p(map(arg), self.c)
        if <PyObject *>x == <PyObject *>None: # c++ pointer compare error
            return Polynomial.__call__(self, *args, **kwds)
        else:
            zz_pX_eval(fx.x, self.x, x.x)
            return self._parent(int(fx))



cdef class Polynomial_dense_modn_ntl_ZZ(Polynomial_dense_mod_n):

    def __init__(self, parent, v=None, check=True, is_gen=False, construct=False):
        if isinstance(v, Polynomial):
            if (<Element>v)._parent == parent:
                Polynomial.__init__(self, parent, is_gen=is_gen)
                self.x = (<Polynomial_dense_modn_ntl_ZZ>v).x
                self.c = (<Polynomial_dense_modn_ntl_ZZ>v).c
                return

        Polynomial_dense_mod_n.__init__(self, parent, v, check=check, is_gen=is_gen, construct=construct)
        cdef ntl_ZZ_pX ntl = self.__poly
        self.__poly = None # this will eventually go away
        self.x = ntl.x
        self.c = ntl.c

    def __dealloc__(self):
        if <object>self.c is not None:
            self.c.restore_c()
        ZZ_pX_destruct(&self.x)

    def ntl_set_directly(self, v):
        # TODO: Get rid of this
        Polynomial_dense_mod_n.ntl_set_directly(self, v)
        # verbatim from __init__
        cdef ntl_ZZ_pX ntl = self.__poly
        self.__poly = None # this will eventually go away
        self.x = ntl.x
        self.c = ntl.c

    cdef Polynomial_dense_modn_ntl_ZZ _new(self):
        cdef Polynomial_dense_modn_ntl_ZZ y = <Polynomial_dense_modn_ntl_ZZ>PY_NEW(Polynomial_dense_modn_ntl_ZZ)
        y.c = self.c
        y._parent = self._parent
        return y

#    def int_list(self):
#        """
#        Returns the coefficents of self as efficiently as possible as a
#        list of python ints.
#
#        EXAMPLES:
#            sage: R.<x> = Integers(100)[]
#            sage: f = x^3 + 5
#            sage: f.int_list()
#            [5, 0, 0, 1]
#            sage: [type(a) for a in f.int_list()]
#            [<type 'long'>, <type 'long'>, <type 'long'>, <type 'long'>]
#        """
#        cdef long i
#        return [ zz_p_rep(zz_pX_GetCoeff(self.x, i)) for i from 0 <= i <= zz_pX_deg(self.x) ]

    def list(self):
        return [self._parent._base(self[n]) for n from 0 <= n <= self.degree()]

    def __getitem__(self, n):
        """
        EXAMPLES:
            sage: R.<x> = Integers(10^30)[]
            sage: f = (x+2)^7
            sage: f[3]
            560
        """
        R = self._parent._base
        if n < 0 or n > ZZ_pX_deg(self.x):
            return R(0)

        self.c.restore_c()
        cdef Integer z
#        cdef ZZ_c rep = ZZ_p_rep(ZZ_pX_coeff(self.x, n))
#        print ZZ_to_int(&rep)
#        ZZ_to_mpz(&z.value, &rep) # does this work?
# TODO, make this faster
        cdef ntl_ZZ_p ntl = ntl_ZZ_p(0, self.c)
        ntl.x = ZZ_pX_coeff(self.x, n)
        return R(ntl._integer_())

    def _unsafe_mutate(self, n, value):
        self.c.restore_c()
        cdef Integer a
        if PY_TYPE_CHECK(value, Integer):
            a = <Integer>value
        else:
            a = ZZ(value)
        cdef ntl_ZZ_p val = ntl_ZZ_p(a, self.c)
        ZZ_pX_SetCoeff(self.x, n, val.x)

    def __getslice__(self, i, j):
        """
        Returns the monomials of self of degree i <= n < j.

        EXAMPLES:
            sage: R.<x> = Integers(10^30)[]
            sage: f = (x+2)^7
            sage: f[3:6]
            84*x^5 + 280*x^4 + 560*x^3
            sage: f[-5:50] == f
            True
            sage: f[6:]
            x^7 + 14*x^6
        """
        R = self.base_ring()
        if i < 0:
            i = 0
        if j > ZZ_pX_deg(self.x)+1:
            j = ZZ_pX_deg(self.x)+1
        v = [ self[t] for t from i <= t < j ]
        return Polynomial_dense_modn_ntl_ZZ(self._parent, v, check=False) << i

    cdef ModuleElement _add_c_impl(self, ModuleElement _right):
        """
        TESTS:
            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) + (x^2 - 6)
            x^2 + x + 999999999999999999999999999999
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e7
        if do_sig: _sig_on
        self.c.restore_c()
        ZZ_pX_add(r.x, self.x, right.x)
        if do_sig: _sig_off
        return r

    cdef ModuleElement _sub_c_impl(self, ModuleElement _right):
        """
        TESTS:
            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) - (x^2 - 6)
            999999999999999999999999999999*x^2 + x + 11
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e7
        if do_sig: _sig_on
        self.c.restore_c()
        ZZ_pX_sub(r.x, self.x, right.x)
        if do_sig: _sig_off
        return r

    cdef RingElement _mul_c_impl(self, RingElement _right):
        """
        TESTS:
            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) * (x^2 - 1)
            x^3 + 5*x^2 + 999999999999999999999999999999*x + 999999999999999999999999999995
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e5
        if do_sig: _sig_on
        self.c.restore_c()
        if self is right:
            ZZ_pX_sqr(r.x, self.x)
        else:
            ZZ_pX_mul(r.x, self.x, right.x)
        if do_sig: _sig_off
        return r

    cdef ModuleElement _rmul_c_impl(self, RingElement c):
        """
        TESTS:
            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) * 3
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = ZZ_pX_deg(self.x) * self.c.p_bits > 1e7
        if do_sig: _sig_on
        self.c.restore_c()
        cdef ntl_ZZ_p value = ntl_ZZ_p(c)
        ZZ_pX_rmul(r.x, self.x, value.x)
        if do_sig: _sig_off
        return r

    cdef ModuleElement _lmul_c_impl(self, RingElement c):
        """
        TESTS:
            sage: R.<x> = Integers(10^30)[]
            sage: 3 * (x+5)
            3*x + 15
        """
        return self._rmul_c_impl(c)

    def __pow__(Polynomial_dense_modn_ntl_ZZ self, ee, dummy):
        """
        TESTS:
            sage: R.<x> = Integers(10^30)[]
            sage: (x+1)^5
            x^5 + 5*x^4 + 10*x^3 + 10*x^2 + 5*x + 1
        """
        cdef bint recip = 0, do_sig
        cdef long e = ee
        if e != ee:
            raise TypeError, "Only integral powers defined."
        elif e < 0:
            recip = 1 # delay because powering frac field elements is slow
            e = -e
        if not self:
            if e == 0:
                raise ArithmeticError, "0^0 is undefined."
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        self.c.restore_c()
        if ZZ_pX_IsX(self.x):
            ZZ_pX_LeftShift(r.x, self.x, e - 1)
        else:
            do_sig = ZZ_pX_deg(self.x) * e * self.c.p_bits > 1e5
            if do_sig: _sig_on
            ZZ_pX_power(r.x, self.x, e)
            if do_sig: _sig_off
        if recip:
            return ~r
        else:
            return r

    def quo_rem(self, right):
        """
        Returns $q$ and $r$, with the degree of $r$ less than the degree of $right$,
        such that $q * right + r = self$.

        EXAMPLES:
            sage: R.<x> = Integers(10^30)[]
            sage: f = x^5+1; g = (x+1)^2
            sage: q, r = f.quo_rem(g)
            sage: q
            x^3 + 999999999999999999999999999998*x^2 + 3*x + 999999999999999999999999999996
            sage: r
            5*x + 5
            sage: q*g + r
            x^5 + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or self._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self.quo_rem(right)
        cdef Polynomial_dense_modn_ntl_ZZ q = self._new()
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(denom.x)) * self.c.p_bits > 1e4
        if do_sig: _sig_on
        self.c.restore_c()
        ZZ_pX_DivRem(q.x, r.x, self.x, denom.x)
        if do_sig: _sig_off
        return q, r

    def __floordiv__(self, right):
        """
        Returns the whole part of self/right, without remainder.

        For q = n // d, we have deg(n - q*d) < deg(d)

        EXAMPLES:
            sage: R.<x> = Integers(10^30)[]
            sage: f = x^7 + 1; g = x^2 - 1
            sage: q = f // g; q
            x^5 + x^3 + x
            sage: f - q*g
            x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self // right
        cdef Polynomial_dense_modn_ntl_ZZ numer = <Polynomial_dense_modn_ntl_ZZ>self
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        cdef Polynomial_dense_modn_ntl_ZZ q = numer._new()
        cdef bint do_sig = (ZZ_pX_deg(numer.x) + ZZ_pX_deg(denom.x)) * numer.c.p_bits > 1e4
        if do_sig: _sig_on
        numer.c.restore_c()
        ZZ_pX_div(q.x, numer.x, denom.x)
        if do_sig: _sig_off
        return q

    def __mod__(self, right):
        """
        EXAMPLES:
            sage: R.<x> = Integers(9^30)[]
            sage: f = x^7 + x + 1; g = x^3 - 1
            sage: r = f % g; r
            2*x + 1
            sage: g * (x^4 + x) + r
            x^7 + x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self % right
        cdef Polynomial_dense_modn_ntl_ZZ numer = <Polynomial_dense_modn_ntl_ZZ>self
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        cdef Polynomial_dense_modn_ntl_ZZ r = numer._new()
        cdef bint do_sig = (ZZ_pX_deg(numer.x) + ZZ_pX_deg(denom.x)) * numer.c.p_bits > 1e4
        if do_sig: _sig_on
        numer.c.restore_c()
        ZZ_pX_rem(r.x, numer.x, denom.x)
        if do_sig: _sig_off
        return r

    def shift(self, n):
        """
        Shift self to left by $n$, which is multiplication by $x^n$,
        truncating if $n$ is negative.

        EXAMPLES:
            sage: R.<x> = Integers(12^30)[]
            sage: f = x^7 + x + 1
            sage: f.shift(1)
            x^8 + x^2 + x
            sage: f.shift(-1)
            x^6 + 1
            sage: f.shift(10).shift(-10) == f
            True
        """
        return self << n

    def __lshift__(Polynomial_dense_modn_ntl_ZZ self, long n):
        """
        TEST:
            sage: R.<x> = Integers(14^30)[]
            sage: f = x^5 + 2*x + 1
            sage: f << 3
            x^8 + 2*x^4 + x^3
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_LeftShift(r.x, self.x, n)
        return r

    def __rshift__(Polynomial_dense_modn_ntl_ZZ self, long n):
        """
        TEST:
            sage: R.<x> = Integers(15^30)[]
            sage: f = x^5 + 2*x + 1
            sage: f >> 3
            x^2
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_RightShift(r.x, self.x, n)
        return r


    def _derivative(self, var=None):
        r"""
        Returns the formal derivative of self with respect to var.

        var must be either the generator of the polynomial ring to which
        this polynomial belongs, or None (either way the behaviour is the
        same).

        SEE ALSO:
            self.derivative()

        EXAMPLES:
            sage: R.<x> = Integers(12^29)[]
            sage: f = x^4 + x + 5
            sage: f._derivative()
            4*x^3 + 1
            sage: f._derivative(None)
            4*x^3 + 1

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

            sage: y = var("y")
            sage: f._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var is not self._parent.gen():
            raise ValueError, "cannot differentiate with respect to %s" % var

        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_diff(r.x, self.x)
        return r


    def reverse(self):
        """
        Reverses the coeffients of self. The reverse of f(x) is x^n f(1/x).

        The degree will go down if the constant term is zero.

        EXAMPLES:
            sage: R.<x> = Integers(12^29)[]
            sage: f = x^4 + 2*x + 5
            sage: f.reverse()
            5*x^4 + 2*x^3 + 1
            sage: f = x^3 + x
            sage: f.reverse()
            x^2 + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_reverse(r.x, self.x)
        return r

    def is_gen(self):
        return ZZ_pX_IsX(self.x)

    def valuation(self):
        """
        Returns the valuation of self, that is, the power of the
        lowest non-zero monomial of self.

        EXAMPLES:
            sage: R.<x> = Integers(10^50)[]
            sage: x.valuation()
            1
            sage: f = x-3; f.valuation()
            0
            sage: f = x^99; f.valuation()
            99
            sage: f = x-x; f.valuation()
            +Infinity
        """
        cdef long n
        cdef ZZ_p_c coeff
        for n from 0 <= n <= ZZ_pX_deg(self.x):
            coeff = ZZ_pX_coeff(self.x, n)
            if not ZZ_p_IsZero(coeff):
                return n
        return infinity

    def __nonzero__(self):
        """
        TESTS:
            sage: R.<x> = Integers(12^29)[]
            sage: f = x^4 + 1
            sage: not f
            False
            sage: not (x-x)
            True
        """
        return not ZZ_pX_IsZero(self.x)

    def degree(self):
        """
        EXAMPLES:
            sage: R.<x> = Integers(14^34)[]
            sage: f = x^4 - x - 1
            sage: f.degree()
            4
            sage: f = 14^43*x + 1
            sage: f.degree()
            0
        """
        return ZZ_pX_deg(self.x)

    def truncate(self, long n):
        """
        Returns this polynomial mod $x^n$.

        EXAMPLES:
            sage: R.<x> = Integers(15^30)[]
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_trunc(r.x, self.x, n)
        return r

    def __call__(self, *args, **kwds):
        """
        Evaluate self at x. If x is a single argument coercable into
        the basering of self, this is done directly in NTL, otherwise
        the generic Polynomial call code is used.

        EXAMPLES:
            sage: R.<x> = Integers(10^30)[]
            sage: f = x^3+7
            sage: f(5)
            132
            sage: f(5r)
            132
            sage: f(mod(5, 10^50))
            132
            sage: f(x)
            x^3 + 7
            sage: S.<y> = Integers(5)[]
            sage: f(y)
            y^3 + 2
        """
        if len(args) != 1 or len(kwds) != 0:
            return Polynomial.__call__(self, *args, **kwds)
        arg = args[0]
        cdef ntl_ZZ_p fx = ntl_ZZ_p(0, self.c), x = None
        if PY_TYPE_CHECK(arg, int) or PY_TYPE_CHECK(arg, Integer):
            x = ntl_ZZ_p(arg, self.c)
        elif PY_TYPE_CHECK(arg, Element):
            if <void *>self._parent._base == <void *>(<Element>arg)._parent: # c++ pointer hack
                x = ntl_ZZ_p(arg, self.c)
            else:
                map = self._parent._base.coerce_map_from((<Element>arg)._parent)
                if map is not None:
                    x = ntl_ZZ_p(map(arg), self.c)
        if <PyObject *>x == <PyObject *>None: # c++ pointer compare error
            return Polynomial.__call__(self, *args, **kwds)
        else:
            ZZ_pX_eval(fx.x, self.x, x.x)
            return self._parent(fx._integer_())


cdef class Polynomial_dense_mod_p(Polynomial_dense_mod_n):
    """
    A dense polynomial over the integers modulo p, where p is prime.
    """

    def gcd(self, right):
        return self._gcd(right)

    def _gcd(self, right):
        """
        Return the GCD of self and other, as a monic polynomial.
        """
        if not isinstance(right, Polynomial_dense_mod_p):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.ntl_ZZ_pX().gcd(right.ntl_ZZ_pX())
        return self.parent()(g, construct=True)

    def xgcd(self, right):
        r"""
        Return the extended gcd of self and other, i.e., elements $r, s, t$ such that
        $$
           r = s \cdot self + t \cdot other.
        $$
        """
        # copied from sage.structure.element.PrincipalIdealDomainElement due to lack of mult inheritance
        if not PY_TYPE_CHECK(right, Element) or not ((<Element>right)._parent is self._parent):
            from sage.rings.arith import xgcd
            return bin_op(self, right, xgcd)
        return self._xgcd(right)

    def _xgcd(self, right):
        """
        Return $g, u, v$ such that \code{g = u*self + v*right}.
        """
        r, s, t = self.ntl_ZZ_pX().xgcd(right.ntl_ZZ_pX())
        return self.parent()(r, construct=True), self.parent()(s, construct=True), \
               self.parent()(t, construct=True)

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial
        OUTPUT:
            an element of the base ring of the polynomial ring

        EXAMPLES:
            sage: R.<x> = GF(19)['x']
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            11
            sage: r.parent() is GF(19)
            True
        """
##        self.parent()._ntl_set_modulus()
        other = self.parent()._coerce_(other)
        return self.base_ring()(str(self.ntl_ZZ_pX().resultant(other.ntl_ZZ_pX())))

    def discriminant(self):
        """
        EXAMPLES:
            sage: _.<x> = PolynomialRing(GF(19))
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            12
        """
##        self.parent()._ntl_set_modulus()
        return self.base_ring()(str(self.ntl_ZZ_pX().discriminant()))

    # PARI is way better than NTL for poly factor for certain degrees, and is called
    # by default in the base class.
    #def factor(self, verbose=False):
    #    M = self.monic()
    #    self.parent()._ntl_set_modulus()
    #    F = [(self.parent()(f, construct=True), n) for f, n in M.ntl_ZZ_pX().factor(verbose)]
    #    return factorization.Factorization(F)
