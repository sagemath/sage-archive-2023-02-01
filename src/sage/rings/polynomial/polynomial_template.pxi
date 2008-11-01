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
from sage.rings.integer cimport Integer
from sage.libs.all import pari_gen

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
        cdef celement gen, monomial, coeff
        cdef Py_ssize_t deg

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            celement_gen(&self.x, 0, parent)

        elif PY_TYPE_CHECK(x, Polynomial_template):
            try:
                celement_set(&self.x, &(<Polynomial_template>x).x, parent)
            except NotImplementedError:
                raise TypeError("%s not understood."%x)

        elif PY_TYPE_CHECK(x, int) or PY_TYPE_CHECK(x, Integer):
            try:
                celement_set_si(&self.x, int(x), parent)
            except NotImplementedError:
                raise TypeError("%s not understood."%x)

        elif PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple):
            parent = (<Polynomial_template>self)._parent

            celement_set_si(&self.x, 0, parent)
            celement_gen(&gen, 0, parent)

            deg = 0
            for e in x:
                # r += parent(e)*power
                celement_pow(&monomial, &gen, deg, NULL, parent)
                coeff = (<Polynomial_template>parent(e)).x
                celement_mul(&monomial, &coeff, &monomial, parent)
                celement_add(&self.x, &self.x, &monomial, parent)
                deg += 1

        elif PY_TYPE_CHECK(x, pari_gen):
            k = (<Polynomial_template>self)._parent.base_ring()
            x = [k(w) for w in x.Vecrev()]
            Polynomial_template.__init__(self, parent, x, check=True, is_gen=False, construct=construct)
        elif PY_TYPE_CHECK(x, Polynomial):
            k = (<Polynomial_template>self)._parent.base_ring()
            x = [k(w) for w in list(x)]
            Polynomial_template.__init__(self, parent, x, check=True, is_gen=False, construct=construct)
        else:
            raise TypeError("Coercion from parent %s not supported."%x.parent())

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
        return [self[i] for i in range(celement_len(&self.x,(<Polynomial_template>self)._parent))]

    def __cinit__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
        """
        celement_construct(&self.x, (<Polynomial_template>self)._parent)

    def __dealloc__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: del x
        """
        celement_destruct(&self.x, (<Polynomial_template>self)._parent)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x + 1
            x + 1
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_add(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._parent)
        return r

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x - 1
            x + 1
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_add(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._parent)
        return r

    def __neg__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: -x
            x
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_neg(&r.x, &self.x, (<Polynomial_template>self)._parent)
        return r

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x*(x+1)
            x^2 + x
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_mul(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._parent)
        return r

    def gcd(self, Polynomial_template other):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x*(x+1)
            sage: f.gcd(x+1)
            x + 1
            sage: f.gcd(x^2)
            x
        """
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_gcd(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>other).x, (<Polynomial_template>self)._parent)
        return r

    def __floordiv__(self, right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x//(x + 1)
            1
            sage: (x + 1)//x
            1
        """
        right = (<Polynomial_template>self)._parent._coerce_(right)
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_floordiv(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._parent)
        return r

    def __mod__(self, other):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: (x^2 + 1) % x^2
            1
        """
        other = (<Polynomial_template>self)._parent._coerce_(other)
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_mod(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>other).x, (<Polynomial_template>self)._parent)
        return r

    def quo_rem(self, right):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x^2 + x + 1
            sage: f.quo_rem(x + 1)
            (x, 1)
        """
        right = (<Polynomial_template>self)._parent._coerce_(right)
        cdef Polynomial_template q = <Polynomial_template>PY_NEW(self.__class__)
        q._parent = (<Polynomial_template>self)._parent
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent
        celement_quorem(&q.x, &r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._parent)
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
        if celement_len(&self.x, (<Polynomial_template>self)._parent) > 1:
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
        return not celement_is_zero(&self.x, (<Polynomial_template>self)._parent)

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
        return celement_cmp(&(<Polynomial_template>left).x, &(<Polynomial_template>right).x, (<Polynomial_template>left)._parent)

    def __hash__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: {x:1}
            {x: 1}
        """
        return celement_hash(&self.x, (<Polynomial_template>self)._parent)

    def __len__(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: P.<x> = GF(2)[]
            sage: len(x)
            2
            sage: len(x+1)
            2
        """
        return celement_len(&self.x, (<Polynomial_template>self)._parent)

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
        cdef long e = ee
        if e != ee:
            raise TypeError("Only integral powers defined.")
        elif e < 0:
            recip = 1 # delay because powering frac field elements is slow
            e = -e
        if not self:
            if e == 0:
                raise ArithmeticError, "0^0 is undefined."
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = (<Polynomial_template>self)._parent

        if modulus is None:
            celement_pow(&r.x, &(<Polynomial_template>self).x, e, NULL, (<Polynomial_template>self)._parent)
        else:
            modulus = (<Polynomial_template>self)._parent._coerce_(modulus)
            celement_pow(&r.x, &(<Polynomial_template>self).x, e, &(<Polynomial_template>modulus).x, (<Polynomial_template>self)._parent)
        if recip:
            return ~r
        else:
            return r

    def is_gen(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x.is_gen()
            True
            sage: (x+1).is_gen()
            False
        """
        cdef celement gen
        celement_gen(&gen, 0, (<Polynomial_template>self)._parent)
        return celement_equal(&self.x, &gen, (<Polynomial_template>self)._parent)

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
        cdef celement gen
        cdef Polynomial_template r
        if n == 0:
            return self

        parent = (<Polynomial_template>self)._parent
        celement_gen(&gen, 0, parent)
        celement_pow(&gen, &gen, abs(n), NULL, parent)
        r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = parent

        if n > 0:
            celement_mul(&r.x, &self.x, &gen, parent)
        else:
            celement_floordiv(&r.x, &self.x, &gen, parent)
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
        cdef celement gen
        cdef Polynomial_template r
        if n == 0:
            return self
        elif n < 0:
            return self >> -n

        parent = (<Polynomial_template>self)._parent
        celement_gen(&gen, 0, parent)
        celement_pow(&gen, &gen, n, NULL, parent)
        r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = parent
        celement_mul(&r.x, &(<Polynomial_template>self).x, &gen, parent)
        return r

    def __rshift__(self, int n):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
        """
        if not PY_TYPE_CHECK(self, Polynomial_template):
            raise TypeError("Cannot %s >> %n."%(self, n))
        cdef celement gen
        cdef Polynomial_template r
        if n == 0:
            return self
        elif n < 0:
            return self >> -n

        parent = (<Polynomial_template>self)._parent
        celement_gen(&gen, 0, parent)
        celement_pow(&gen, &gen, n, NULL, parent)
        r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = parent

        celement_floordiv(&r.x, &(<Polynomial_template>self).x, &gen, parent)
        return r

    def is_zero(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
        """
        return celement_is_zero(&self.x, (<Polynomial_template>self)._parent)

    def is_one(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
        """
        return celement_is_one(&self.x, (<Polynomial_template>self)._parent)

    def degree(self):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
        """
        return Integer(celement_len(&self.x, (<Polynomial_template>self)._parent)-1)

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
        if n < 0:
            raise ValueError(" n must be >= 0.")
        parent = (<Polynomial_template>self)._parent
        cdef celement gen
        celement_gen(&gen, 0, parent)
        celement_pow(&gen, &gen, n, NULL, parent)

        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = parent
        celement_mod(&r.x, &self.x, &gen, parent)
        return r
