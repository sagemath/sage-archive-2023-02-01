"""
Polynomial Template for C/C++ Library Interfaces

AUTHOR:
    -- Robert Bradshaw (2008-10) original idea for templating
    -- Martin Albrecht (2008-10) initial implementation

This file implements a simple templating engine for linking univariate
polynomials to their C/C++ library implementations. It requires a
'linkage' file which implements the \code{celement_} functions (see
\code{sage.libs.ntl.ntl_GF2X_linkage} for an example). Both parts are
then pluygged together by inclusion of the linkage file when
inheriting from this class. See
\code{sage.rings.polynomial.polynomial_gf2x} for an example.

We illustrate the generic glueing using univariate polynomials over
GF(2).

NOTE:
Implementations using this template MUST implement coercion from base
ring elements and \code{__getitem__}. See \code{Polynomial_GF2X} for an example.
"""
#*****************************************************************************
#       Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#       Copyright (C) 2008 Robert Bradshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.element cimport ModuleElement, Element, RingElement
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.integer cimport Integer
from sage.libs.all import pari_gen

from sage.interfaces.all import singular as singular_default

def make_element(parent, args):
    return parent(*args)

cdef class Polynomial_template(Polynomial):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: P(0)
            0
            sage: P(GF(2)(1))
            1
            sage: P(3)
            1
            sage: P([1,0,1])
            x^2 + 1
            sage: P(map(GF(2),[1,0,1]))
            x^2 + 1
        """
        cdef celement *gen, *monomial
        cdef Py_ssize_t deg
        cdef cparent _parent

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            celement_construct(&self.x, get_cparent(parent))
            celement_gen(&self.x, 0, get_cparent(parent))

        elif PY_TYPE_CHECK(x, Polynomial_template):
            try:
                celement_construct(&self.x, get_cparent(parent))
                celement_set(&self.x, &(<Polynomial_template>x).x, get_cparent(parent))
            except NotImplementedError:
                raise TypeError("%s not understood."%x)

        elif PY_TYPE_CHECK(x, int) or PY_TYPE_CHECK(x, Integer):
            try:
                celement_construct(&self.x, get_cparent(parent))
                celement_set_si(&self.x, int(x), get_cparent(parent))
            except NotImplementedError:
                raise TypeError("%s not understood."%x)

        elif PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple):
            parent = (<Polynomial_template>self)._parent
            _parent = get_cparent(parent)

            celement_construct(&self.x, get_cparent(parent))
            gen = celement_new(_parent)
            monomial = celement_new(_parent)

            celement_set_si(&self.x, 0, _parent)
            celement_gen(gen, 0, _parent)

            deg = 0
            for e in x:
                # r += parent(e)*power
                celement_pow(monomial, gen, deg, NULL, _parent)
                celement_mul(monomial, &(<Polynomial_template>self.__class__(parent, e)).x, monomial, _parent)
                celement_add(&self.x, &self.x, monomial, _parent)
                deg += 1

            celement_delete(gen, _parent)
            celement_delete(monomial, _parent)

        elif PY_TYPE_CHECK(x, dict):
            parent = (<Polynomial_template>self)._parent
            _parent = get_cparent(parent)

            celement_construct(&self.x, get_cparent(parent))

            gen = celement_new(_parent)
            monomial = celement_new(_parent)

            celement_set_si(&self.x, 0, _parent)
            celement_gen(gen, 0, _parent)

            for deg, coef in x.iteritems():
                celement_pow(monomial, gen, deg, NULL, _parent)
                celement_mul(monomial, &(<Polynomial_template>parent(coef)).x, monomial, _parent)
                celement_add(&self.x, &self.x, monomial, _parent)

            celement_delete(gen, _parent)
            celement_delete(monomial, _parent)

        elif PY_TYPE_CHECK(x, pari_gen):
            k = (<Polynomial_template>self)._parent.base_ring()
            x = [k(w) for w in x.Vecrev()]
            self.__class__.__init__(self, parent, x, check=True, is_gen=False, construct=construct)
        elif PY_TYPE_CHECK(x, Polynomial):
            k = (<Polynomial_template>self)._parent.base_ring()
            x = [k(w) for w in list(x)]
            Polynomial_template.__init__(self, parent, x, check=True, is_gen=False, construct=construct)
        elif PY_TYPE_CHECK(x, FractionFieldElement) and (x.parent().base() is parent or x.parent().base() == parent) and x.denominator() == 1:
            x = x.numerator()
            self.__class__.__init__(self, parent, x, check=check, is_gen=is_gen, construct=construct)
        else:
            x = parent.base_ring()(x)
            self.__class__.__init__(self, parent, x, check=check, is_gen=is_gen, construct=construct)

    def __reduce__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: loads(dumps(x)) == x
            True
        """
        return make_element, ((<Polynomial_template>self)._parent, (self.list(), False, self.is_gen()))

    def list(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x.list()
            [0, 1]
            sage: list(x)
            [0, 1]
        """
        cdef Py_ssize_t i
        cdef cparent _parent = get_cparent((<Polynomial_template>self)._parent)
        return [self[i] for i in range(celement_len(&self.x, _parent))]

    def __dealloc__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: del x
        """
        celement_destruct(&self.x, get_cparent((<Polynomial_template>self)._parent))

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x + 1
            x + 1
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        celement_add(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, get_cparent((<Polynomial_template>self)._parent))
        #assert(r._parent(pari(self) + pari(right)) == r)
        return r

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x - 1
            x + 1
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        celement_sub(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, get_cparent((<Polynomial_template>self)._parent))
        #assert(r._parent(pari(self) - pari(right)) == r)
        return r

    def __neg__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: -x
            x
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        celement_neg(&r.x, &self.x, get_cparent((<Polynomial_template>self)._parent))
        #assert(r._parent(-pari(self)) == r)
        return r

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x*(x+1)
            x^2 + x
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        celement_mul(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, get_cparent((<Polynomial_template>self)._parent))
        #assert(r._parent(pari(self) * pari(right)) == r)
        return r

    def gcd(self, Polynomial_template other):
        """
        Return the greatest common divisor of self and other.

        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x*(x+1)
            sage: f.gcd(x+1)
            x + 1
            sage: f.gcd(x^2)
            x
        """
        if(celement_is_zero(&self.x, get_cparent((<Polynomial_template>self)._parent))):
            return other
        if(celement_is_zero(&other.x, get_cparent((<Polynomial_template>self)._parent))):
            return self

        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        celement_gcd(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>other).x, get_cparent((<Polynomial_template>self)._parent))
        #assert(r._parent(pari(self).gcd(pari(other))) == r)
        return r

    def xgcd(self, Polynomial_template other):
        """
        Computes extended gcd of self and other.

        EXAMPLE:
            sage: P.<x> = GF(7)[]
            sage: f = x*(x+1)
            sage: f.xgcd(x+1)
            (x + 1, 0, 1)
            sage: f.xgcd(x^2)
            (x, 1, 6)
        """
        if(celement_is_zero(&self.x, get_cparent((<Polynomial_template>self)._parent))):
            return other, self._parent(0), self._parent(1)
        if(celement_is_zero(&other.x, get_cparent((<Polynomial_template>self)._parent))):
            return self, self._parent(1), self._parent(0)

        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent

        cdef Polynomial_template s = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&s.x, get_cparent((<Polynomial_template>self)._parent))
        s._parent = (<Polynomial_template>self)._parent

        cdef Polynomial_template t = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&t.x, get_cparent((<Polynomial_template>self)._parent))
        t._parent = (<Polynomial_template>self)._parent

        celement_xgcd(&r.x, &s.x, &t.x, &(<Polynomial_template>self).x, &(<Polynomial_template>other).x, get_cparent((<Polynomial_template>self)._parent))
        #rp, sp, tp = pari(self).xgcd(pari(other))
        #assert(r._parent(rp) == r)
        #assert(s._parent(sp) == s)
        #assert(t._parent(tp) == t)
        return r,s,t

    def __floordiv__(self, right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x//(x + 1)
            1
            sage: (x + 1)//x
            1
        """
        cdef Polynomial_template _right = <Polynomial_template>(<Polynomial_template>self)._parent._coerce_(right)
        if celement_is_zero(&_right.x, get_cparent((<Polynomial_template>self)._parent)):
            raise ZeroDivisionError
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        #assert(r._parent(pari(self) // pari(right)) == r)
        celement_floordiv(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, get_cparent((<Polynomial_template>self)._parent))
        return r

    def __mod__(self, other):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: (x^2 + 1) % x^2
            1
        """
        cdef Polynomial_template _other = <Polynomial_template>(<Polynomial_template>self)._parent._coerce_(other)
        if celement_is_zero(&_other.x, get_cparent((<Polynomial_template>self)._parent)):
            raise ZeroDivisionError

        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        celement_mod(&r.x, &(<Polynomial_template>self).x, &_other.x, get_cparent((<Polynomial_template>self)._parent))
        #assert(r._parent(pari(self) % pari(other)) == r)
        return r

    def quo_rem(self, right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x^2 + x + 1
            sage: f.quo_rem(x + 1)
            (x, 1)
        """
        cdef Polynomial_template _right = <Polynomial_template>(<Polynomial_template>self)._parent._coerce_(right)

        if celement_is_zero(&_right.x, get_cparent((<Polynomial_template>self)._parent)):
            raise ZeroDivisionError

        cdef Polynomial_template q = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&q.x, get_cparent((<Polynomial_template>self)._parent))
        q._parent = (<Polynomial_template>self)._parent

        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent

        celement_quorem(&q.x, &r.x, &(<Polynomial_template>self).x, &_right.x, get_cparent((<Polynomial_template>self)._parent))
        return q,r

    def __long__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: int(x)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int

            sage: int(P(1))
            1
        """
        if celement_len(&self.x, get_cparent((<Polynomial_template>self)._parent)) > 1:
            raise ValueError("Cannot coerce polynomial with degree %d to integer."%(self.degree()))
        return int(self[0])

    def __nonzero__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: bool(x), x.is_zero()
            (True, False)
            sage: bool(P(0)), P(0).is_zero()
            (False, True)
        """
        return not celement_is_zero(&self.x, get_cparent((<Polynomial_template>self)._parent))

    def __richcmp__(left, right, int op):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x != 1
            True
            sage: x < 1
            False
            sage: x > 1
            True
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
        """
        return celement_cmp(&(<Polynomial_template>left).x, &(<Polynomial_template>right).x, get_cparent((<Polynomial_template>left)._parent))

    def __hash__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: {x:1}
            {x: 1}
        """
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        cdef long c_hash
        cdef long var_name_hash
        cdef int i
        for i from 0<= i <= self.degree():
            if i == 1:
                # we delay the hashing until now to not waste it one a constant poly
                var_name_hash = hash(self.variable_name())
            #  I'm assuming (incorrectly) that hashes of zero indicate that the element is 0.
            # This assumption is not true, but I think it is true enough for the purposes and it
            # it allows us to write fast code that omits terms with 0 coefficients.  This is
            # important if we want to maintain the '==' relationship with sparse polys.
            c_hash = hash(self[i])
            if c_hash != 0:
                if i == 0:
                    result += c_hash
                else:
                    # Hash (self[i], generator, i) as a tuple according to the algorithm.
                    result_mon = c_hash
                    result_mon = (1000003 * result_mon) ^ var_name_hash
                    result_mon = (1000003 * result_mon) ^ i
                    result += result_mon
        if result == -1:
            return -2
        return result


    def __pow__(self, ee, modulus):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: P.<x> = GF(2)[]
            sage: x^1000
            x^1000
            sage: (x+1)^2
            x^2 + 1
            sage: (x+1)^(-2)
            1/(x^2 + 1)
            sage: f = x^9 + x^7 + x^6 + x^5 + x^4 + x^2 + x
            sage: h = x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + 1
            sage: (f^2) % h
            x^9 + x^8 + x^7 + x^5 + x^3
            sage: pow(f, 2, h)
            x^9 + x^8 + x^7 + x^5 + x^3
        """
        if not PY_TYPE_CHECK(self, Polynomial_template):
            raise NotImplementedError("%s^%s not defined."%(ee,self))
        cdef bint recip = 0, do_sig

        cdef long e
        try:
            e = ee
        except OverflowError:
            return Polynomial.__pow__(self, ee, modulus)
        if e != ee:
            raise TypeError("Only integral powers defined.")
        elif e < 0:
            recip = 1 # delay because powering frac field elements is slow
            e = -e
        if not self:
            if e == 0:
                raise ArithmeticError, "0^0 is undefined."
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent

        if modulus is None:
            celement_pow(&r.x, &(<Polynomial_template>self).x, e, NULL, get_cparent((<Polynomial_template>self)._parent))
        else:
            modulus = (<Polynomial_template>self)._parent._coerce_(modulus)
            celement_pow(&r.x, &(<Polynomial_template>self).x, e, &(<Polynomial_template>modulus).x, get_cparent((<Polynomial_template>self)._parent))

        #assert(r._parent(pari(self)**ee) == r)
        if recip:
            return ~r
        else:
            return r

    def __copy__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: copy(x) is x
            False
            sage: copy(x) == x
            True
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, get_cparent((<Polynomial_template>self)._parent))
        r._parent = (<Polynomial_template>self)._parent
        celement_set(&r.x, &self.x, get_cparent((<Polynomial_template>self)._parent))
        return r

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
        cdef cparent _parent = get_cparent((<Polynomial_template>self)._parent)
        if i < 0:
            i = 0
        if j > celement_len(&self.x,_parent):
            j = celement_len(&self.x, _parent)
        x = (<Polynomial_template>self)._parent.gen()
        v = [self[t] for t from i <= t < j]

        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        Polynomial_template.__init__(r, (<Polynomial_template>self)._parent, v)
        return r << i

    def is_gen(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x.is_gen()
            True
            sage: (x+1).is_gen()
            False
        """
        cdef celement *gen = celement_new(get_cparent((<Polynomial_template>self)._parent))
        celement_gen(gen, 0, get_cparent((<Polynomial_template>self)._parent))
        cdef bint r = celement_equal(&self.x, gen, get_cparent((<Polynomial_template>self)._parent))
        celement_delete(gen, get_cparent((<Polynomial_template>self)._parent))
        return r

    def shift(self, int n):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1
            sage: f.shift(1)
            x^4 + x^3 + x
            sage: f.shift(-1)
            x^2 + x
        """
        cdef Polynomial_template r
        if n == 0:
            return self

        parent = (<Polynomial_template>self)._parent
        cdef cparent _parent = get_cparent(parent)
        cdef celement *gen = celement_new(_parent)
        celement_gen(gen, 0, _parent)
        celement_pow(gen, gen, abs(n), NULL, _parent)
        r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, _parent)
        r._parent = parent

        if n > 0:
            celement_mul(&r.x, &self.x, gen, _parent)
        else:
            celement_floordiv(&r.x, &self.x, gen, _parent)
        celement_delete(gen, _parent)
        return r

    def __lshift__(self, int n):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1
            sage: f << 1
            x^4 + x^3 + x
        """
        if not PY_TYPE_CHECK(self, Polynomial_template):
            raise TypeError("Cannot %s << %n."%(self, n))
        cdef celement *gen, *tmp
        cdef Polynomial_template r
        if n == 0:
            return self
        elif n < 0:
            return self >> -n

        parent = (<Polynomial_template>self)._parent
        cdef cparent _parent = get_cparent(parent)
        gen = celement_new(_parent)
        celement_gen(gen, 0, _parent)
        celement_pow(gen, gen, n, NULL, _parent)
        r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, _parent)
        r._parent = parent
        celement_mul(&r.x, &(<Polynomial_template>self).x, gen, _parent)
        celement_delete(gen, _parent)
        return r

    def __rshift__(self, int n):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x>>1
            1
            sage: (x^2 + x)>>1
            x + 1
            sage: (x^2 + x) >> -1
            x^3 + x^2
        """
        if not PY_TYPE_CHECK(self, Polynomial_template):
            raise TypeError("Cannot %s >> %n."%(self, n))
        cdef celement *gen
        cdef Polynomial_template r
        if n == 0:
            return self
        elif n < 0:
            return self << -n

        parent = (<Polynomial_template>self)._parent
        cdef cparent _parent = get_cparent(parent)
        gen = celement_new(_parent)
        celement_gen(gen, 0, _parent)
        celement_pow(gen, gen, n, NULL, _parent)
        r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, _parent)
        r._parent = parent

        celement_floordiv(&r.x, &(<Polynomial_template>self).x, gen, _parent)
        celement_delete(gen, _parent)
        return r

    def is_zero(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x.is_zero()
            False
        """
        return celement_is_zero(&self.x, get_cparent((<Polynomial_template>self)._parent))

    def is_one(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: P(1).is_one()
            True
        """
        return celement_is_one(&self.x, get_cparent((<Polynomial_template>self)._parent))

    def degree(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x.degree()
            1
            sage: P(1).degree()
            0
            sage: P(0).degree()
            -1
        """
        return Integer(celement_len(&self.x, get_cparent((<Polynomial_template>self)._parent))-1)

    cpdef Polynomial truncate(self, long n):
        """
        Returns this polynomial mod $x^n$.

        EXAMPLES:
            sage: R.<x> =GF(2)[]
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        parent = (<Polynomial_template>self)._parent
        cdef cparent _parent = get_cparent(parent)

        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        celement_construct(&r.x, _parent)
        r._parent = parent

        if n <= 0:
            return r

        cdef celement *gen = celement_new(_parent)
        celement_gen(gen, 0, _parent)
        celement_pow(gen, gen, n, NULL, _parent)

        celement_mod(&r.x, &self.x, gen, _parent)
        celement_delete(gen, _parent)
        return r

    def _singular_(self, singular=singular_default, have_ring=False, force=False):
        r"""
        Return \Singular representation of this polynomial

        INPUT:
            singular -- \Singular interpreter (default: default interpreter)
            have_ring -- set to True if the ring was already set in \Singular
            force -- ignored.

        EXAMPLE:
            sage: P.<x> = PolynomialRing(GF(7))
            sage: f = 3*x^2 + 2*x + 5
            sage: singular(f)
            3*x^2+2*x-2
        """
        if not have_ring:
            self.parent()._singular_(singular,force=force).set_ring() #this is expensive
        return singular(self._singular_init_())

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

        P = self.parent()
        x = P.gen()
        res = P(0)
        for i,c in enumerate(self.list()[1:]):
            res += (i+1)*c*x**i
        return res
