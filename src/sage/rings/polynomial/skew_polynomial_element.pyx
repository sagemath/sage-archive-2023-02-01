r"""
This module implements elements in skew polynomial rings.

DEFINITION::

Let `R` be a commutative ring equipped with an endomorphism `\sigma`.

The skew polynomial ring over `(R, \sigma)` is the ring `S = `R[X,\sigma]`
is the usual ring of polynomials over `R` equipped with the skew
multiplication defined by the rule `X*a = \sigma(a)*X` for all `a`
in `R`.

EXAMPLES::

We illustrate some functionalities implemented in this class.

We create the skew polynomial ring::

    sage: R.<t> = ZZ[]
    sage: sigma = R.hom([t+1])
    sage: S.<x> = R['x',sigma]; S
    Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

and some elements in it::

    sage: a = t + x + 1; a
    x + t + 1
    sage: b = S([t^2,t+1,1]); b
    x^2 + (t + 1)*x + t^2
    sage: c = S.random_element(degree=3,monic=True); c
    x^3 + (-95*t^2 + t + 2)*x^2 + (-t^2 + t)*x + 2*t - 8

Ring operations are supported::

    sage: a + b
    x^2 + (t + 2)*x + t^2 + t + 1
    sage: a - b
    -x^2 - t*x - t^2 + t + 1

    sage: a * b
    x^3 + (2*t + 3)*x^2 + (2*t^2 + 4*t + 2)*x + t^3 + t^2
    sage: b * a
    x^3 + (2*t + 4)*x^2 + (2*t^2 + 3*t + 2)*x + t^3 + t^2
    sage: a * b == b * a
    False

    sage: b^2
    x^4 + (2*t + 4)*x^3 + (3*t^2 + 7*t + 6)*x^2 + (2*t^3 + 4*t^2 + 3*t + 1)*x + t^4
    sage: b^2 == b*b
    True


Sage implements also some arithmetics over skew polynomial rings.
You will find below a short panorama.

DEFINITION:

Let `a` and `b` be two skew polynomials in the same ring `S`.
The *left (resp. right) euclidean division* of `a` by `b` is
a couple `(q,r)` of elements in `S` such that

-  `a = q*b + r` (resp. `a = b*q + r`)

-  the degree of `r` is less than the degree of `b`

`q` (resp. `r`) is called the *quotient* (resp. the remainder)
of this euclidean division.

PROPERTIES:

Keeping the previous notations, if the leading coefficient of `b`
is a unit (e.g. if `b` is monic) then the quotient and the remainder
in the *right* euclidean division exist and are unique.

The same result holds for the *left* euclidean division if in addition
the twist map defining the skew polynomial ring is invertible.

EXAMPLES::

    sage: q,r = c.quo_rem(b)   # default side is right
    sage: q
    x - 95*t^2
    sage: r
    (95*t^3 + 93*t^2 - t - 1)*x + 95*t^4 + 2*t - 8
    sage: c == q*b + r
    True

The operators ``//`` and ``%`` give respectively the quotient
and the remainder of the right euclidean division::

    sage: q == c // b
    True
    sage: r == c % b
    True

If we want left euclidean division, we need to specify
``side=Left``. Nonetheless, it won't work over our current
`S` because Sage can't invert the twist map::

    sage: q,r = c.quo_rem(b,side=Left)
    Traceback (most recent call last):
    ...
    NotImplementedError

Here is a working example over a finite field::

    sage: k.<t> = GF(5^3)
    sage: Frob = k.frobenius_endomorphism()
    sage: S.<x> = k['x',Frob]
    sage: a = x^4 + (4*t + 1)*x^3 + (t^2 + 3*t + 3)*x^2 + (3*t^2 + 2*t + 2)*x + 3*t^2 + 3*t + 1
    sage: b = (2*t^2 + 3)*x^2 + (3*t^2 + 1)*x + 4*t + 2
    sage: q,r = a.quo_rem(b,side=Left)
    sage: q
    (4*t^2 + t + 1)*x^2 + (2*t^2 + 2*t + 2)*x + 2*t^2 + 4*t + 3
    sage: r
    (t + 2)*x + 3*t^2 + 2*t + 4
    sage: a == b*q + r
    True


Once we have euclidean divisions, we have for free gcd and lcm
(at least if the base ring is a field).
This class provides an implementation of gcd and lcm::

    sage: a = (x + t) * (x + t^2)^2
    sage: b = (x + t) * (t*x + t + 1) * (x + t^2)
    sage: a.gcd(b)  # default side is right
    x + t^2
    sage: a.gcd(b,side=Left)
    x + t

For lcm, the default side is left but be very careful: by
convention, a left (resp. right) lcm is common multiple on
the right (resp. left)::

    sage: c = a.lcm(b); c  # default side is left
    x^5 + (4*t^2 + t + 3)*x^4 + (3*t^2 + 4*t)*x^3 + 2*t^2*x^2 + (2*t^2 + t)*x + 4*t^2 + 4
    sage: c.is_divisible_by(a)
    True
    sage: c.is_divisible_by(b)
    True

    sage: d = a.lcm(b,side=Right); d
    x^5 + (t^2 + 1)*x^4 + (3*t^2 + 3*t + 3)*x^3 + (3*t^2 + t + 2)*x^2 + (4*t^2 + 3*t)*x + 4*t + 4
    sage: d.is_divisible_by(a,side=Left)
    True
    sage: d.is_divisible_by(b,side=Left)
    True

AUTHOR:

-  Xavier Caruso (2012-06-29)
"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************


include "../../ext/stdsage.pxi"

import operator, copy, re

import skew_polynomial_ring
import sage.rings.infinity as infinity
from sage.misc.latex import latex
from sage.structure.factorization import Factorization

from sage.categories.homset import Hom

from sage.structure.element import RingElement
from sage.structure.element cimport Element, RingElement, ModuleElement

from sage.rings.ring import Field

from sage.structure.parent_gens cimport ParentWithGens

from sage.rings.integer cimport Integer
from sage.categories.map cimport Map
from sage.rings.morphism cimport RingHomomorphism
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

from sage.structure.side import Left, Right

def is_SkewPolynomial(a):
    """
    Return True if `a` is a skew polynomial (over some base).

    INPUT:

    -  ``a`` -- an object

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_element import is_SkewPolynomial
        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]
        sage: a = x^2
        sage: is_SkewPolynomial(a)
        True
    """
    return type(a) is SkewPolynomial
#    return PY_TYPE_CHECK(a, SkewPolynomial)


cdef class CenterSkewPolynomial_generic_dense(Polynomial_generic_dense):
    """
    A class for elements in the center of a skew polynomial ring.
    """
    pass


cdef class SkewPolynomial(AlgebraElement):
    """
    A skew polynomial.
    """
    def __init__(self,parent,is_gen=False,construct=False):
        """
        The following examples illustrate the creation of elements of
        skew polynomial rings.

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: P = x+t; P
            x + t
            sage: Q = S([1,t,t+2]); Q
            (t + 2)*x^2 + t*x + 1
        """
        AlgebraElement.__init__(self, parent)
        self._is_gen = is_gen

    # you may have to replicate this boilerplate code in derived classes if you override
    # __richcmp__.  The python documentation at  http://docs.python.org/api/type-structs.html
    # explains how __richcmp__, __hash__, and __cmp__ are tied together.
    def __hash__(self):
        return self._hash_c()

    cdef long _hash_c(self):
        """
        This hash incorporates the name of the variable.
        """
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        cdef long c_hash
        cdef long var_name_hash
        cdef int i
        for i from 0<= i <= self.degree():
            if i == 1:
                # we delay the hashing until now to not waste it one a constant poly
                var_name_hash = hash((<ParentWithGens>self._parent)._names[0])
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


    # Comparison
    # ----------

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp_(right, op)

    cpdef int _cmp_(self, Element other) except -2:
        """
        Compare the two skew polynomials self and other.

        We order polynomials first by degree, then in dictionary order
        starting with the coefficient of largest degree.
        """
        cdef list x = (<SkewPolynomial>self)._list_c()
        cdef list y = (<SkewPolynomial>other)._list_c()
        cdef Py_ssize_t dx = len(x), dy = len(y)
        cdef Py_ssize_t i
        cdef int c
        c = cmp(dx,dy)
        if c: return c
        for i from dx > i >= 0:
            c = cmp(x[i],y[i])
            if c: return c
        return 0


    # Some c functions
    # ----------------

    cdef SkewPolynomial _new_c(self,list coeffs,Parent P,char check=0):
        """
        Fast creation of a new skew polynomial

        .. NOTE::

            Override this function in classes which inherit
            from SkewPolynomial.
        """
        return P(list)


    cpdef SkewPolynomial _new_constant_poly(self,RingElement a,Parent P,char check=0):
        """
        Fast creation of a new constant skew polynomial
        """
        if a:
            return self._new_c([a],P,check)
        else:
            return self._new_c([],P)


    cdef list _list_c(self):
        """
        Return the list of the underlying elements of this
        skew polynomial.

        .. WARNING::

            It is *not* a copy; do not modify this list!
        """
        raise NotImplementedError


    # A skew polynomial as a list of coefficients
    # -------------------------------------------

    def list(self):
        """
        Return a new copy of the list of the underlying elements of self.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: l = a.list(); l
            [t^2 + 1, 0, t + 1, 0, 1]

        Note that v is a list, it is mutable, and each call to the list
        method returns a new list::

            sage: type(l)
            <type 'list'>
            sage: l[0] = 5
            sage: a.list()
            [t^2 + 1, 0, t + 1, 0, 1]
        """
        return list((<SkewPolynomial>self)._list_c())


    def __getitem__(self,n):
        """
        INPUT::

        -  ``n`` -- an integer

        OUTPUT::

        -  the `n`-th coefficient of self

        .. NOTE::

            Coefficients are on the left (ie before the variable)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^2 + (t + 3/7)*x + t^2
            sage: a[1]
            t + 3/7
            sage: a[3]
            0
        """
        try:
            return (<SkewPolynomial>self)._list_c()[n]
        except IndexError:
            return self.base_ring()(0)


    def __getslice__(self, Py_ssize_t i, Py_ssize_t j):
        """
        EXAMPLES::
            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^2 + (t + 3/7)*x + t^2
            sage: a[1:]
            t*x^2 + (t + 3/7)*x
            sage: a[:1]
            t^2
            sage: a[3:]
            0
        """
        if i <= 0:
            i = 0
            zeros = []
        elif i > 0:
            zeros = [self._parent.base_ring()(0)] * i
        return self._new_c(zeros + self._list_c()[i:j], self._parent, 1)


    def __setitem__(self, n, value):
        """
        Set the `n`-th coefficient of this skew polynomial. This always
        raises an IndexError, since in Sage polynomials are immutable.

        INPUT:

        -  ``n`` - an integer

        -  ``value`` - value to set the `n`-th coefficient to

        OUTPUT: an IndexError is always raised.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: a[1] = t + 1
            Traceback (most recent call last):
            ...
            IndexError: skew polynomials are immutable
        """
        raise IndexError, "skew polynomials are immutable"


    # Basic operations
    # ----------------

    def degree(self):
        """
        Return the degree of this skew polynomial. The zero
        skew polynomial has degree -1.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t*x^3 + t^2*x + 1
            sage: a.degree()
            3

        By convention, the degree of 0 is -1::

            sage: S(0).degree()
            -1
        """
        return len((<SkewPolynomial>self)._list_c())-1


    def valuation(self):
        """
        Return the valuation of this skew polynomial.
        The zero skew polynomial has valuation +Infinity.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t*x^3 + t^2*x
            sage: a.valuation()
            1

        By convention, the valuation of 0 is +Infinity::

            sage: S(0).valuation()
            +Infinity
        """
        cdef list x = (<SkewPolynomial>self)._list_c()
        if len(x) == 0:
            return infinity.infinity
        cdef Py_ssize_t v = 0
        while x[v].is_zero():
            v += 1
        return v


    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Compute self + right

        INPUT:

        -  right -- a skew polynomial over the same base

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(monic=True); a
            x^2 + (-12*t^2 + 1/2*t - 1/95)*x - 1/2*t^2 - 4
            sage: b = -S.random_element(monic=True); b
            -x^2 + (5/2*t - 2/3)*x + 1/4*t^2 - 1/2*t + 1
            sage: c = a+b; c
            (-12*t^2 + 3*t - 193/285)*x - 1/4*t^2 - 1/2*t - 3
            sage: c.degree()
            1
        """
        cdef Py_ssize_t i, min
        cdef list x = (<SkewPolynomial>self)._list_c()
        cdef list y = (<SkewPolynomial>right)._list_c()
        cdef Py_ssize_t dx = len(x), dy = len(y)

        if dx > dy:
            return self._new_c([x[i] + y[i] for i from 0 <= i < dy] + x[dy:], self._parent, 0)
        elif dx < dy:
            return self._new_c([x[i] + y[i] for i from 0 <= i < dx] + y[dx:], self._parent, 0)
        else:
            return self._new_c([x[i] + y[i] for i from 0 <= i < dx], self._parent, 1)


    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Compute self - right

        INPUT:

        -  right -- a skew polynomial over the same base

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(monic=True); a
            x^2 + (-12*t^2 + 1/2*t - 1/95)*x - 1/2*t^2 - 4
            sage: b = S.random_element(monic=True); b
            x^2 + (-5/2*t + 2/3)*x - 1/4*t^2 + 1/2*t - 1
            sage: c = a-b; c
            (-12*t^2 + 3*t - 193/285)*x - 1/4*t^2 - 1/2*t - 3
            sage: c.degree()
            1
        """
        cdef Py_ssize_t i, min
        cdef list x = (<SkewPolynomial>self)._list_c()
        cdef list y = (<SkewPolynomial>right)._list_c()
        cdef Py_ssize_t dx = len(x), dy = len(y)
        cdef RingElement c

        if dx > dy:
            return self._new_c([x[i] - y[i] for i from 0 <= i < dy] + x[dy:], self._parent, 0)
        elif dx < dy:
            return self._new_c([x[i] - y[i] for i from 0 <= i < dx] + [ -c for c in y[dx:] ], self._parent, 0)
        else:
            return self._new_c([x[i] - y[i] for i from 0 <= i < dx], self._parent, 1)


    cpdef ModuleElement _neg_(self):
        """
        Return the opposite of self

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^2 + x - 3
            sage: -a
            -t*x^2 - x + 3
        """
        return self._new_c([-x for x in (<SkewPolynomial>self)._list_c()], self._parent, 0)


    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Compute self * right (in this order)

        INPUT:

        -  right -- an element of the base ring

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t
            sage: b = t
            sage: a * b
            (t + 1)*x + t^2
            sage: a * b == b * a
            False
        """
        if right == 0:
            return self._parent(0)
        cdef list x = (<SkewPolynomial>self)._list_c()
        cdef Py_ssize_t i
        map = self._parent._map
        return self._new_c([ (map**i)(right)*x[i] for i from 0 <= i < len(x) ], self._parent, 0)


    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        Compute left * self (in this order)

        INPUT:

        -  left -- an element of the base ring

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t
            sage: b = x + t
            sage: a * b
            t*x + t^2
            sage: a * b == b * a
            False
        """
        if left == 0:
            return self.parent().zero()
#            return self.parent().zero_element()
        cdef list x = (<SkewPolynomial>self)._list_c()
        cdef Py_ssize_t i
        map = self._parent._map
        return self._new_c([ left*x[i] for i from 0 <= i < len(x) ], self._parent, 0)


    cpdef RingElement _mul_(self, RingElement right):
        """
        Compute self * right (in this order)

        INPUT:

        -  right -- a skew polynomial in the same ring

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t; a
            x^2 + t
            sage: b = x^2 + (t + 1)*x; b
            x^2 + (t + 1)*x
            sage: a * b
            x^4 + (t + 3)*x^3 + t*x^2 + (t^2 + t)*x
            sage: a * b == b * a
            False
        """
        cdef list x = (<SkewPolynomial>self)._list_c()
        cdef list y = (<SkewPolynomial>right)._list_c()
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t dx = len(x)-1, dy = len(y)-1
        parent = self._parent
        if dx == -1:
            return self
        elif dy == -1:
            return right
        elif dx == 0:
            c = x[0]
            return self._new_c([c*a for a in y], parent, 0)
        cdef list coeffs = []
        for k from 0 <= k <= dx+dy:
            start = 0 if k <= dy else k-dy # max(0, k-dy)
            end = k if k <= dx else dx # min(k, dx)
            sum = x[start] * parent.twist_map(start)(y[k-start])
            for i from start < i <= end:
                sum += x[i] * parent.twist_map(i)(y[k-i])
            coeffs.append(sum)
        return self._new_c(coeffs, parent, 0)


    def square(self):
        """
        Return the square of self

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t; a
            x + t
            sage: a.square()
            x^2 + (2*t + 1)*x + t^2
            sage: a.square() == a*a
            True
        """
        return self * self


    def conjugate(self,n):
        """
        INPUT:

        -  `n` -- an integer

        OUTPUT:

        -  this skew polynomial conjugated by x^n (where x is
           the variable)

        .. NOTE::

            This conjugate is obtained from the skew polynomial by
            applying the `n`-th iterate of the twist map to each of
            its coefficients.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^3 + (t^2 + 1)*x^2 + 2*t
            sage: b = a.conjugate(2); b
            (t + 2)*x^3 + (t^2 + 4*t + 5)*x^2 + 2*t + 4
            sage: x^2*a == b*x^2
            True

        In principle, negative values for `n` are allowed... but Sage
        needs to be able to invert the twist map::

            sage: b = a.conjugate(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        Here is a working example::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: u = T.random_element(); u
            (2*t^2 + 3)*y^2 + (4*t^2 + t + 4)*y + 2*t^2 + 2
            sage: v = u.conjugate(-1); v
            (3*t^2 + t)*y^2 + (4*t^2 + 2*t + 4)*y + 3*t^2 + t + 4
            sage: u*y == y*v
            True
        """
        return self._new_c([ self._parent.twist_map(n)(x) for x in (<SkewPolynomial>self)._list_c() ], self._parent, 0)


    # Other useful mathematical functions
    # -----------------------------------

    def constant_coefficient(self):
        """
        Return the constant coefficient (i.e. the coefficient of degree
        `0`) of this skew polynomial.

        OUTPUT: element of base ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t^2 + 2
            sage: a.constant_coefficient()
            t^2 + 2
        """
        cdef x = (<SkewPolynomial>self)._list_c()
        if len(x) == 0:
            return self.base_ring()(0)
        else:
            return x[0]


    def leading_coefficient(self):
        """
        Return the leading coefficient of this skew polynomial.

        OUTPUT: element of the base ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + (t+1)*x^5 + t^2*x^3 - x^5
            sage: a.leading_coefficient()
            t
        """
        cdef x = (<SkewPolynomial>self)._list_c()
        if len(x) == 0:
            raise ValueError("")
        return x[-1]


    def __nonzero__(self):
        return self.degree() >= 0


    def is_unit(self):
        """
        Return True if this skew polynomial is a unit.

        Not yet implemented.
        """
        raise NotImplementedError


    def is_monic(self):
        """
        Returns True if this skew polynomial is monic. The zero polynomial
        is by definition not monic.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t
            sage: a.is_monic()
            True
            sage: a = 0*x
            sage: a.is_monic()
            False
            sage: a = t*x^3 + x^4 + (t+1)*x^2
            sage: a.is_monic()
            True
            sage: a = (t^2 + 2*t)*x^2 + x^3 + t^10*x^5
            sage: a.is_monic()
            False
        """
        return not self.is_zero() and self[self.degree()] == 1


    def lmonic(self):
        """
        Return the unique monic skew polynomial `a` of the same
        degree which divides this skew polynomial on the left

        .. Note::

            This skew polynomial is self dividing on the
            *right* by the `n`-th iterative (`n` is the degree of
            self) of the inverse of the twist map applied to the
            leading coefficient.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = a.lmonic(); b
            x^3 + (4*t^2 + 3*t)*x^2 + (4*t + 2)*x + 2*t^2 + 4*t + 3

        Check list::

            sage: b.degree() == a.degree()
            True
            sage: b.is_divisible_by(a,side=Left)
            True
            sage: twist = S.twist_map(-a.degree())
            sage: a == b * twist(a.leading_coefficient())
            True

        Note that `b` does not divise `a` on the right::

            sage: b.is_divisible_by(a,side=Right)
            False

        This function does not work if the leading coefficient is not a
        unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x
            sage: a.lmonic()
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient is not a unit
        """
        try:
            a = self.base_ring()(~self.leading_coefficient())
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient is not a unit")
        return self*self._parent.twist_map(-self.degree())(a)


    def rmonic(self):
        """
        Return the unique monic skew polynomial `a` of the same
        degree which divides this skew polynomial on the right

        .. Note::

            This skew polynomial is self dividing on the *left*
            by its leading coefficient.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = a.rmonic(); b
            x^3 + (2*t^2 + 3*t + 4)*x^2 + (3*t^2 + 4*t + 1)*x + 2*t^2 + 4*t + 3

        Check list::

            sage: b.degree() == a.degree()
            True
            sage: b.is_divisible_by(a,side=Right)
            True
            sage: a == a.leading_coefficient() * b
            True

        Note that `b` does not divise `a` on the right::

            sage: b.is_divisible_by(a,side=Left)
            False

        This function does not work if the leading coefficient is not a
        unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x
            sage: a.rmonic()
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient is not a unit
        """
        try:
            a = self.base_ring()(~self.leading_coefficient())
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient is not a unit")
        return a*self


    def monic(self, side=Right):
        """
        INPUT:

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        -  the unique monic skew polynomial `a` of the same
           degree which divides this skew polynomial on ``side``

        .. Note::

            if ``side`` is ``Right``, this skew polynomial is self
            dividing on the *left* by its leading coefficient; if
            ``side`` is ``Left``, it is self dividing on the
            *right* by the `n`-th iterative (`n` is the degree of
            self) of the inverse of the twist map applied to the
            leading coefficient.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = a.monic(); b
            x^3 + (2*t^2 + 3*t + 4)*x^2 + (3*t^2 + 4*t + 1)*x + 2*t^2 + 4*t + 3

        Check list::

            sage: b.degree() == a.degree()
            True
            sage: b.is_divisible_by(a,side=Right)
            True
            sage: a == a.leading_coefficient() * b
            True

        Note that `b` does not divise `a` on the left::

            sage: b.is_divisible_by(a,side=Left)
            False

        The same on the left::

            sage: b = a.monic(side=Left); b
            x^3 + (4*t^2 + 3*t)*x^2 + (4*t + 2)*x + 2*t^2 + 4*t + 3
            sage: b.degree() == a.degree()
            True
            sage: b.is_divisible_by(a,side=Left)
            True
            sage: twist = S.twist_map(-a.degree())
            sage: a == b * twist(a.leading_coefficient())
            True

        This function does not work if the leading coefficient is not a
        unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x
            sage: a.monic()
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient is not a unit
        """
        if self.is_monic():
            return self
        if side == Right:
            return self.rmonic()
        else:
            return self.lmonic()


    # Divisibility
    # ------------


    def lquo_rem(self, other):
        """
        INPUT:

        -  ``other`` -- a skew polynomial ring over the same
           base ring

        OUTPUT:

        -  the quotient and the remainder of the left euclidean
           division of this skew polynomial by ``other``

        .. NOTE::

            Doesn't work if the leading coefficient of the divisor
            is not a unit or if Sage can't invert the twist map.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = (3*t^2 + 4*t + 2)*x^2 + (2*t^2 + 4*t + 3)*x + 2*t^2 + t + 1
            sage: q,r = a.lquo_rem(b)
            sage: a == b*q + r
            True

        In the following example, Sage does not know the inverse
        of the twist map::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (-2*t^2 - t + 1)*x^3 + (-t^2 + t)*x^2 + (-12*t - 2)*x - t^2 - 95*t + 1
            sage: b = x^2 + (5*t - 6)*x - 4*t^2 + 4*t - 1
            sage: a.lquo_rem(b)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        cdef list a = list((<SkewPolynomial>self)._list_c())
        cdef list b = (<SkewPolynomial>other)._list_c()
        cdef Py_ssize_t i, j
        cdef Py_ssize_t da = self.degree(), db = other.degree()
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            return self._new_c([],self._parent), self
        try:
            inv = self.base_ring()(~b[db])
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient of the divisor is not invertible")
        cdef list q = [ ]
        parent = self._parent
        for i from da-db >= i >= 0:
            c = parent.twist_map(-db)(inv*a[i+db])
            for j from 0 <= j < db:
                a[i+j] -= b[j] * parent.twist_map(j)(c)
            q.append(c)
        q.reverse()
        return self._new_c(q,parent), self._new_c(a[:db],parent,1)


    def rquo_rem(self, other):
        """
        INPUT:

        -  ``other`` -- a skew polynomial ring over the same
           base ring

        OUTPUT:

        -  the quotient and the remainder of the left euclidean
           division of this skew polynomial by ``other``

        .. NOTE::

            Doesn't work if the leading coefficient of the divisor
            is not a unit.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(degree=4); a
            t^2*x^4 + (-12*t^2 - 2*t - 1)*x^3 + (-95*t^2 + t + 2)*x^2 + (-t^2 + t)*x + 2*t - 8
            sage: b = S.random_element(monic=True); b
            x^2 + (4*t^2 - t - 2)*x - t^2 + t - 1
            sage: q,r = a.rquo_rem(b)
            sage: a == q*b + r
            True

        The leading coefficient of the divisor need to be invertible::

            sage: c = S.random_element(); c
            (-4*t^2 + t)*x^2 - 2*t^2*x + 5*t^2 - 6*t - 4
            sage: a.rquo_rem(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        cdef list a = list((<SkewPolynomial>self)._list_c())
        cdef list b = (<SkewPolynomial>other)._list_c()
        cdef Py_ssize_t i, j
        cdef Py_ssize_t da = self.degree(), db = other.degree()
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            return self._new_c([],parent), self
        try:
            inv = self.base_ring()(~b[db])
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient of the divisor is not invertible")
        cdef list q = [ ]
        parent = self._parent
        for i from da-db >= i >= 0:
            c = parent.twist_map(i)(inv) * a[i+db]
            for j from 0 <= j < db:
                a[i+j] -= c * parent.twist_map(i)(b[j])
            q.append(c)
        q.reverse()
        return self._new_c(q,parent), self._new_c(a[:db],parent,1)


    def quo_rem(self,other,side=Right):
        r"""
        INPUT:

        -  ``other`` -- a skew polynomial ring over the same
           base ring

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        -  the quotient and the remainder of the euclidean
           division of this skew polynomial by ``other``

        .. NOTE::

            Doesn't work if the leading coefficient of the divisor
            is not a unit.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(degree=4); a
            t^2*x^4 + (-12*t^2 - 2*t - 1)*x^3 + (-95*t^2 + t + 2)*x^2 + (-t^2 + t)*x + 2*t - 8
            sage: b = S.random_element(monic=True); b
            x^2 + (4*t^2 - t - 2)*x - t^2 + t - 1
            sage: q,r = a.quo_rem(b)  # right euclidean division
            sage: a == q*b + r
            True

        The left euclidean division doesn't work over this `S` because
        Sage cannot invert `\sigma`::

            sage: q,r = a.quo_rem(b,side=Left)
            Traceback (most recent call last):
            ...
            NotImplementedError

        In any case, the leading coefficient of the divisor need to be
        invertible::

            sage: c = S.random_element(); c
            (-4*t^2 + t)*x^2 - 2*t^2*x + 5*t^2 - 6*t - 4
            sage: a.quo_rem(c)  # right euclidean division
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible

        When the base ring is a finite field, everything works fine::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = 4*x^3 + (t^2 + 2*t + 4)*x^2 + (4*t^2 + 2*t + 4)*x + t^2 + 2
            sage: b = (3*t + 2)*x^2 + 2*t*x + 3*t^2
            sage: q,r = a.quo_rem(b)  # right euclidean division
            sage: a == q*b + r
            True
            sage: q,r = a.quo_rem(b,side=Left)
            sage: a == b*q + r
            True
        """
        if side == Right:
            return self.rquo_rem(other)
        else:
            return self.lquo_rem(other)


    def rem(self,other,side=Right):
        r"""
        INPUT:

        -  ``other`` -- a skew polynomial ring over the same
           base ring

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        -  the remainder of the euclidean division of this
           skew polynomial by ``other``

        .. NOTE::

            Doesn't work if the leading coefficient of the divisor
            is not a unit.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: b = x^2 + 2*t*x + 2
            sage: a = (x+t)*b + t*x + 1
            sage: a.rem(b)   # right euclidean division
            t*x + 1

        The left euclidean division doesn't work over this `S` because
        Sage cannot invert `\sigma`::

            sage: a.rem(b,side=Left)
            Traceback (most recent call last):
            ...
            NotImplementedError

        In any case, the leading coefficient of the divisor need to be
        invertible::

            sage: (a*t).rem(b*t)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible

        When the base ring is a finite field, everything works fine::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: b = x^2 + 2*t*x + 2
            sage: a = (x+t)*b + t*x + 1
            sage: a.rem(b)  # right euclidean division
            t*x + 1
            sage: a.rem(b,side=Left)
            (3*t^2 + 2)*x + 2*t^2
        """
        _,r = self.quo_rem(other,side=side)
        return r


    def __mod__(self,other):
        """
        Return the remainder in the *right* euclidean division of
        this skew polynomial by ``other``

        TESTS::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: b = x^2 + 2*t*x + 2
            sage: a = (x+t)*b + t*x + 1
            sage: a % b
            t*x + 1

            sage: (a*t).rem(b*t)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        _,r = self.rquo_rem(other)
        return r


    def __floordiv__(self,right):
        """
        Return the quotient of the right euclidean division of self by right.

        The algorithm fails if the leading coefficient of the divisor (right)
        is not invertible.

        .. SEEALSO::


        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: b = x^2 + t
            sage: a = (x^2 + t*x + 1)*b + t^3*x
            sage: a // b
            x^2 + t*x + 1

            sage: (t*a) // (t*b)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible

        """
        q,_ = self.rquo_rem(right)
        return q


    cpdef RingElement _div_(self, RingElement right):
        """
        Not Implemented (since localization of Ore rings is
        not yet implemented, see trac #13215).

        Use the operator `//` even for exact division.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^5 + (t + 2)*x^2 + t^2
            sage: b = x^3 + 4*t
            sage: c = a*b

            sage: c / b
            Traceback (most recent call last):
            ...
            NotImplementedError: Please use `//` (even for exact division)

            sage: c // b == a
            True
        """
        raise NotImplementedError("Please use `//` (even for exact division)")


    def is_divisible_by(self,other,side=Right):
        """
        INPUT:

        -  ``other`` -- a skew polynomial over the same base

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        Return True iff self is divisible by other on ``side``

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a*b
            sage: c.is_divisible_by(a)   # on the right
            False
            sage: c.is_divisible_by(b)
            True
            sage: c.is_divisible_by(a,side=Left)
            True
            sage: c.is_divisible_by(b,side=Left)
            False

        Divisibility by 0 has no sense::

            sage: c.is_divisible_by(S(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        This function does not work if the leading coefficient of the divisor
        is not a unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + 2*x + t
            sage: b = (t+1)*x + t^2
            sage: c = a*b
            sage: c.is_divisible_by(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        return self.rem(other,side=side).is_zero()


    def divides(self,other,side=Right):
        """
        INPUT:

        -  ``other`` -- a skew polynomial over the same base

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        Return True iff self divides other on ``side``

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a*b
            sage: a.divides(c)   # on the right
            False
            sage: b.divides(c)
            True
            sage: a.divides(c,side=Left)
            True
            sage: b.divides(c,side=Left)
            False

        Divisibility by 0 has no sense::

            sage: S(0).divides(c)
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        This function does not work if the leading coefficient of the divisor
        is not a unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + 2*x + t
            sage: b = (t+1)*x + t^2
            sage: c = a*b
            sage: b.divides(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        return other.rem(self,side=side).is_zero()


    # greastest commun divisor
    # ------------------------

    def lxgcd(self,other,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The left gcd of self and other, that is a skew polynomial
        `g` with the following property: any skew polynomial is
        divisible on the left by `g` iff it is divisible on the left
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        - Two skew polynomials `u` and `v` such that:

        .. MATH::

            g = self * u + other * v

        .. NOTE::

            Works only if two following conditions are fulfilled
            (otherwise left gcd do not exist in general):
            1) the base ring is a field and 2) the twist map on
            this field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: g,u,v = a.lxgcd(b); g
            x + t
            sage: a*u + b*v == g
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: g,u,v = a.lxgcd(b,monic=False); g
            2*t*x + 4*t + 2
            sage: a*u + b*v == g
            True

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.lxgcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map need to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.lxgcd(b)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(self.base_ring(),Field):
            raise TypeError("the base ring must be a field")
        G = self
        U = self._parent(1)
        if other.is_zero():
            V = self._parent(0)
        else:
            V1 = self._parent(0)
            V3 = other
            while not V3.is_zero():
                Q,R = G.lquo_rem(V3)
                T = U - V1*Q
                U = V1
                G = V3
                V1 = T
                V3 = R
            V,_ = (G-self*U).lquo_rem(other)
        if monic:
            lc = ~G.leading_coefficient()
            lc = self._parent.twist_map(-G.degree())(lc)
            G = G*lc
            U = U*lc
            V = V*lc
        return G,U,V


    def rxgcd(self,other,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The right gcd of self and other, that is a skew polynomial
        `g` with the following property: any skew polynomial is
        divisible on the right by `g` iff it is divisible on the right
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        - Two skew polynomials `u` and `v` such that:

        .. MATH::

            g = u * self + v * other

        .. NOTE::

            Works only if the base ring is a field (otherwise right
            gcd do not exist in general).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: g,u,v = a.rxgcd(b); g
            x + t
            sage: u*a + v*b == g
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: g,u,v = a.rxgcd(b,monic=False); g
            (4*t^2 + 4*t + 1)*x + 4*t^2 + 4*t + 3
            sage: u*a + v*b == g
            True

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.rxgcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
        """
        if not isinstance(self.base_ring(),Field):
            raise TypeError("the base ring must be a field")
        G = self
        U = self._parent(1)
        if other.is_zero():
            V = self._parent(0)
        else:
            V1 = self._parent(0)
            V3 = other
            while not V3.is_zero():
                Q, R = G.rquo_rem(V3)
                T = U - Q*V1
                U = V1
                G = V3
                V1 = T
                V3 = R
            V,_ = (G-U*self).rquo_rem(other)
        if monic:
            lc = ~G.leading_coefficient()
            G = lc*G
            U = lc*U
            V = lc*V
        return G,U,V


    def xgcd(self,other,side=Right,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The left (resp. right) gcd of self and other, that is a
        skew polynomial `g` with the following property: any skew
        polynomial is divisible on the left (resp. right) by `g`
        iff it is divisible on the left (resp. right) by both
        ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        - Two skew polynomials `u` and `v` such that:

            - (if side=Left) `g = self * u + other * v`

            - (if side=Right) `g = u * self + v * other`

        .. NOTE::

            Works only if the base ring is a field and, when
            ``side`` is ``Left`` if the twist map on this field
            is bijective (otherwise gcd do not exist in general).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + 2*t) * (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x + 2*t) * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: g,u,v = a.xgcd(b); g    # side = Right
            x + t
            sage: u*a + v*b == g
            True
            sage: g,u,v = a.xgcd(b,side=Left); g
            x + 2*t
            sage: a*u + b*v == g
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: g,u,v = a.xgcd(b,monic=False); g
            (2*t^2 + 2*t + 4)*x + 2*t^2 + 3*t + 4
            sage: u*a + v*b == g
            True
            sage: g,u,v = a.xgcd(b,side=Left,monic=False); g
            (2*t^2 + 3*t + 4)*x + 2*t^2 + t + 3
            sage: a*u + b*v == g
            True

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + 2*t) * (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x + 2*t) * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.xgcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
            sage: a.xgcd(b,side=Left)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map need to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + 2*t) * (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x + 2*t) * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: g,u,v = a.xgcd(b); g
            x + t
            sage: g,u,v = a.xgcd(b,side=Left); g
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if side == Right:
            return self.rxgcd(other,monic=monic)
        else:
            return self.lxgcd(other,monic=monic)


    def rgcd(self,other,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The right gcd of self and other, that is a skew polynomial
        `g` with the following property: any skew polynomial is
        divisible on the right by `g` iff it is divisible on the right
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field (otherwise right
            gcd do not exist in general).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.rgcd(b)
            x + t

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.rgcd(b,monic=False)
            (4*t^2 + 4*t + 1)*x + 4*t^2 + 4*t + 3

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.rgcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
        """
        if not isinstance(self.base_ring(),Field):
            raise TypeError("the base ring must be a field")
        if other.is_zero():
            return self
        A = self
        B = other
        while not B.is_zero():
            A,B = B, A % B
        if monic:
            A = A.rmonic()
        return A


    def lgcd(self,other,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The left gcd of self and other, that is a skew polynomial
        `g` with the following property: any skew polynomial is
        divisible on the left by `g` iff it is divisible on the left
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if two following conditions are fulfilled
            (otherwise left gcd do not exist in general):
            1) the base ring is a field and 2) the twist map on
            this field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.lgcd(b)
            x + t

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.lgcd(b,monic=False)
            2*t*x + 4*t + 2

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.lgcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map need to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.lgcd(b)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(self.base_ring(),Field):
            raise TypeError("the base ring must be a field")
        if other.is_zero():
            return self
        A = self
        B = other
        while not B.is_zero():
            A,B = B, A.rem(B,side=Left)
        if monic:
            A = A.lmonic()
        return A


    def gcd(self,other,side=Right,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The left (resp. right) gcd of self and other, that is a
        skew polynomial `g` with the following property: any skew
        polynomial is divisible on the left (resp. right) by `g`
        iff it is divisible on the left (resp. right) by both
        ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field and, when
            ``side=Left`` if, in addition, the twist map on this
            field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + 2*t) * (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x + 2*t) * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.gcd(b);    # side = Right
            x + t
            sage: a.gcd(b,side=Left)
            x + 2*t

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.gcd(b,monic=False)
            (2*t^2 + 2*t + 4)*x + 2*t^2 + 3*t + 4
            sage: a.gcd(b,side=Left,monic=False)
            (2*t^2 + 3*t + 4)*x + 2*t^2 + t + 3

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + 2*t) * (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x + 2*t) * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.gcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
            sage: a.gcd(b,side=Left)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And, when ``side=Left``, the twist map need to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + 2*t) * (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x + 2*t) * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.gcd(b)
            x + t
            sage: a.gcd(b,side=Left)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if side == Right:
            return self.rgcd(other,monic=monic)
        else:
            return self.lgcd(other,monic=monic)


    # Lowest common multiple
    # ----------------------

    def llcm(self,other,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The left lcm of self and other, that is a skew polynomial
        `g` with the following property: any skew polynomial divides
        `g` on the *right* iff it divides both ``self`` and ``other``
        on the *right*.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field (otherwise left
            lcm do not exist in general).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t^2) * (x + t)
            sage: b = 2 * (x^2 + t + 1) * (x * t)
            sage: c = a.llcm(b); c
            x^5 + (2*t^2 + t + 4)*x^4 + (3*t^2 + 4)*x^3 + (3*t^2 + 3*t + 2)*x^2 + (t^2 + t + 2)*x
            sage: c.is_divisible_by(a,side=Right)
            True
            sage: c.is_divisible_by(b,side=Right)
            True
            sage: a.degree() + b.degree() == c.degree() + a.rgcd(b).degree()
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.llcm(b,monic=False)
            (t^2 + t)*x^5 + (4*t^2 + 4*t + 1)*x^4 + (t + 1)*x^3 + (t^2 + 2)*x^2 + (3*t + 4)*x

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t^2) * (x + t)
            sage: b = 2 * (x^2 + t + 1) * (x * t)
            sage: a.llcm(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
        """
        if not isinstance(self.base_ring(),Field):
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError
        U = self._parent(1)
        G = self
        V1 = self._parent(0)
        V3 = other
        while not V3.is_zero():
            Q, R = G.rquo_rem(V3)
            T = U - Q*V1
            U = V1
            G = V3
            V1 = T
            V3 = R
        V1 = V1*self
        if monic:
            V1 = V1.rmonic()
        return V1


    def rlcm(self,other,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The right lcm of self and other, that is a skew polynomial
        `g` with the following property: any skew polynomial divides
        `g` on the *left* iff it divides both ``self`` and ``other``
        on the *left*.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if two following conditions are fulfilled
            (otherwise right lcm do not exist in general):
            1) the base ring is a field and 2) the twist map on
            this field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t) * (x + t^2)
            sage: b = 2 * (x + t) * (x^2 + t + 1)
            sage: c = a.rlcm(b); c
            x^4 + (2*t^2 + t + 2)*x^3 + (3*t^2 + 4*t + 1)*x^2 + (3*t^2 + 4*t + 1)*x + t^2 + 4
            sage: c.is_divisible_by(a,side=Left)
            True
            sage: c.is_divisible_by(b,side=Left)
            True
            sage: a.degree() + b.degree() == c.degree() + a.lgcd(b).degree()
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.rlcm(b,monic=False)
            2*t*x^4 + (3*t + 1)*x^3 + (4*t^2 + 4*t + 3)*x^2 + (3*t^2 + 4*t + 2)*x + 3*t^2 + 2*t + 3

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t) * (x + t^2)
            sage: b = 2 * (x + t) * (x^2 + t + 1)
            sage: a.rlcm(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map need to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x + t^2)
            sage: b = 2 * (x + t) * (x^2 + t + 1)
            sage: a.rlcm(b)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(self.base_ring(),Field):
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError
        R = self.parent()
#        U = R.one_element()
        U = R.one()
        G = self
#        V1 = R.zero_element()
        V1 = R.zero()
        V3 = other
        while not V3.is_zero():
            Q, R = G.lquo_rem(V3)
            T = U - V1*Q
            U = V1
            G = V3
            V1 = T
            V3 = R
        V1 = self*V1
        if monic:
            V1 = V1.lmonic()
        return V1


    def lcm(self,other,side=Left,monic=True):
        """
        INPUT:

        -  ``other`` -- an other skew polynomial over the same
           base

        -  ``side`` -- ``Left`` or ``Right`` (default: Left)

        -  ``monic`` -- boolean (default: True)

        OUTPUT:

        - The left (resp. right) lcm of self and other, that is a
        skew polynomial `g` with the following property: any skew
        polynomial divides `g` on the right (resp. left) iff it
        divides both ``self`` and ``other`` on the rignt (resp.
        left)
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field and, when
            ``side=Right`` if, in addition, the twist map on this
            field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + 2*t) * (x + t^2) * (x + t)
            sage: b = 2 * (x + 2*t) * (x + t + 1) * (x + t)

            sage: c = a.lcm(b); c      # side = Left
            x^5 + (2*t^2 + 4*t)*x^4 + (2*t^2 + 3*t + 4)*x^3 + (t^2 + t + 3)*x^2 + (2*t^2 + 4*t + 3)*x + t^2 + t + 1
            sage: c.is_divisible_by(a,side=Right)
            True
            sage: c.is_divisible_by(b,side=Right)
            True

            sage: c = a.lcm(b,side=Right); c
            x^5 + (3*t^2 + 2*t + 1)*x^4 + (4*t^2 + t + 2)*x^3 + (4*t^2 + 2*t + 3)*x^2 + 4*t^2*x + t^2 + 2
            sage: c.is_divisible_by(a,side=Left)
            True
            sage: c.is_divisible_by(b,side=Left)
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.lcm(b,monic=False)
            (3*t^2 + 3*t + 4)*x^5 + (2*t^2 + 4*t + 1)*x^4 + (t^2 + t + 1)*x^3 + (2*t^2 + t + 4)*x^2 + (t^2 + 3*t + 3)*x + t^2 + 1
            sage: a.lcm(b,side=Right,monic=False)
            (4*t + 4)*x^5 + (3*t^2 + t + 1)*x^4 + (t^2 + 4)*x^3 + (4*t^2 + 2*t + 4)*x^2 + (3*t^2 + t)*x + 2*t^2 + 2

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + 2*t) * (x + t^2) * (x + t)
            sage: b = 2 * (x + 2*t) * (x + t + 1) * (x + t)
            sage: a.lcm(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
            sage: a.lcm(b,side=Right)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map need to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + 2*t) * (x + t^2) * (x + t)
            sage: b = 2 * (x + 2*t) * (x + t + 1) * (x + t)
            sage: a.lcm(b,side=Right)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if side == Right:
            return self.rlcm(other,monic=monic)
        else:
            return self.llcm(other,monic=monic)


    # Printing
    # --------

    def _repr_(self):
        """
        Return string representation of this skew polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t^2 + 1/2*x*t;
            sage: a._repr_()
            '(1/2*t + 1/2)*x + t^2'
        """
        return self._repr()


    def _repr(self,name=None):
        """
        Return string representation of this skew polynomial.

        INPUT:

        -  ``name`` -- the name of the variable (default: the
           name given when the skew polynomial ring was created)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t^2 + 1/2*x*t;
            sage: a._repr()
            '(1/2*t + 1/2)*x + t^2'
            sage: a._repr(name='y')
            '(1/2*t + 1/2)*y + t^2'
        """
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = self.list()
        for n in reversed(xrange(m)):
            x = coeffs[n]
            if x:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(r' 1(\.0+)?\*',' ', s)
        s = re.sub(r' -1(\.0+)?\*',' -', s)
        if s == " ":
            return "0"
        return s[1:]


    def _latex_(self,name=None):
        """
        Return a latex representation of this skew polynomial.

        INPUT:

        -  ``name`` -- the name of the variable (default: the
           name given when the skew polynomial ring was created)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t^2 + 1/2*x*t;
            sage: a._latex_()
            '\\left(\\frac{1}{2} t + \\frac{1}{2}\\right) x + t^{2}'
            sage: a._latex_(name='y')
            '\\left(\\frac{1}{2} t + \\frac{1}{2}\\right) y + t^{2}'
        """
        s = " "
        coeffs = self.list()
        m = len(coeffs)
        if name is None:
            name = self.parent().latex_variable_names()[0]
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        for n in reversed(xrange(m)):
            x = coeffs[n]
            x = y = latex(x)
            if x != '0':
                if n != m-1:
                    s += " + "
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if n > 1:
                    var = "|%s^{%s}"%(name,n)
                elif n==1:
                    var = "|%s"%name
                else:
                    var = ""
                s += "%s %s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(" 1(\.0+)? \|"," ", s)
        s = re.sub(" -1(\.0+)? \|", " -", s)
        s = s.replace("|","")
        if s==" ":
            return "0"
        return s[1:].lstrip().rstrip()


    # Misc
    # ----

    def _is_atomic(self):
        return (self.degree() == self.valuation() and
                self.leading_coefficient()._is_atomic())


    def base_ring(self):
        """
        Return the base ring of this skew polynomial.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element()
            sage: a.base_ring()
            Univariate Polynomial Ring in t over Integer Ring
            sage: a.base_ring() is R
            True
        """
        return self.parent().base_ring()


    def shift(self, n):
        """
        Return this skew polynomial multiplied on the right by the
        power `x^n`. If `n` is negative, terms below `x^n` will be
        discarded.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^5 + t^4*x^4 + t^2*x^2 + t^10
            sage: a.shift(0)
            x^5 + t^4*x^4 + t^2*x^2 + t^10
            sage: a.shift(-1)
            x^4 + t^4*x^3 + t^2*x
            sage: a.shift(-5)
            1
            sage: a.shift(2)
            x^7 + t^4*x^6 + t^2*x^4 + t^10*x^2

        One can also use the infix shift operator::

            sage: a >> 2
            x^3 + t^4*x^2 + t^2
            sage: a << 2
            x^7 + t^4*x^6 + t^2*x^4 + t^10*x^2
        """
        if n == 0 or self.degree() < 0:
            return self
        if n > 0:
            return self._parent(n*[self.base_ring().zero()] + self.list(), check=False)
#            return self._parent(n*[self.base_ring().zero_element()] + self.list(), check=False)
        if n < 0:
            if n > self.degree():
                return self._parent([])
            else:
                return self._parent(self.list()[-n:], check=False)


    def __lshift__(self, k):
        return self.shift(k)


    def __rshift__(self, k):
        return self.shift(-k)


    def change_variable_name(self, var):
        """
        Return a new polynomial over the same base ring but in a different
        variable.

        INPUT:

        -  ``var`` -- the name of the new variable

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^3 + (2*t + 1)*x  + t^2 + 3*t + 5
            sage: b = a.change_variable_name('y'); b
            y^3 + (2*t + 1)*y  + t^2 + 3*t + 5

        Remark that a new parent is created at the same time::

            sage: b.parent()
            Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        parent = self._parent
        R = parent.base_ring()[var,parent.twist_map()]
        return R(self.list())


    def __copy__(self):
        """
        Return a copy of this skew polynomial.
        """
        return self


    def dict(self):
        """
        Return a sparse dictionary representation of this skew
        polynomial.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2012 + t*x^1006 + t^3 + 2*t
            sage: a.dict()
            {0: t^3 + 2*t, 2012: 1, 1006: t}
        """
        X = {}
        Y = self.list()
        for i in xrange(len(Y)):
            c = Y[i]
            if c:
                X[i] = c
        return X


    def is_monomial(self):
        """
        Returns True if self is a monomial, i.e., a power of the generator.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: x.is_monomial()
            True
            sage: (x+1).is_monomial()
            False
            sage: (x^2).is_monomial()
            True
            sage: S(1).is_monomial()
            True

        The coefficient must be 1::

            sage: (2*x^5).is_monomial()
            False
            sage: S(t).is_monomial()
            False

        To allow a non-1 leading coefficient, use is_term()::

            sage: (2*x^5).is_term()
            True
            sage: S(t).is_term()
            True
        """
        return len(self.exponents()) == 1 and self.leading_coefficient() == 1


    def is_term(self):
        """
        Return True if self is an element of the base ring times a
        power of the generator.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: x.is_term()
            True
            sage: R(1).is_term()
            True
            sage: (3*x^5).is_term()
            True
            sage: (1+3*x^5).is_term()
            False

        To require that the coefficient is 1, use is_monomial() instead::

            sage: (3*x^5).is_monomial()
            False
        """
        return len(self.exponents()) == 1


    def is_gen(self):
        return self._is_gen


    def coefficients(self):
        """
        Return the nonzero coefficients of the monomials
        appearing in self.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: a.coefficients()
            [t^2 + 1, t + 1, 1]
        """
        return [c for c in self.list() if not c.is_zero()]


    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: a.exponents()
            [0, 2, 4]
        """
        l = self.list()
        return [i for i in range(len(l)) if not l[i].is_zero()]


    def prec(self):
        """
        Return the precision of this polynomial. This is always infinity,
        since polynomials are of infinite precision by definition (there is
        no big-oh).

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: x.prec()
            +Infinity
        """
        return infinity.infinity


    def padded_list(self,n=None):
        """
        Return list of coefficients of self up to (but not including)
        `q^n`.

        Includes 0's in the list on the right so that the list has length
        `n`.

        INPUT:


        -  ``n`` - (default: None); if given, an integer that
           is at least 0


        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + t*x^3 + t^2*x^5
            sage: a.padded_list()
            [1, 0, 0, t, 0, t^2]
            sage: a.padded_list(10)
            [1, 0, 0, t, 0, t^2, 0, 0, 0, 0]
            sage: len(a.padded_list(10))
            10
            sage: a.padded_list(3)
            [1, 0, 0]
            sage: a.padded_list(0)
            []
            sage: a.padded_list(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be at least 0
        """
        v = self.list()
        if n is None:
            return v
        if n < 0:
            raise ValueError("n must be at least 0")
        if len(v) < n:
            z = self._parent.base_ring().zero()
#            z = self._parent.base_ring().zero_element()
            return v + [z]*(n - len(v))
        else:
            return v[:int(n)]


    def coeffs(self):
        r"""
        Returns ``self.list()``.

        (It is potentially slightly faster to use
        ``self.list()`` directly.)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: a.coeffs()
            [t^2 + 1, 0, t + 1, 0, 1]
        """
        return self.list()


    def variable_name(self):
        """
        Return name of variable used in this polynomial as a string.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t
            sage: a.variable_name()
            'x'
        """
        return self.parent().variable_name()


cdef class SkewPolynomial_generic_dense(SkewPolynomial):
    def __init__(self, parent, x=None, int check=1, is_gen=False, int construct=0, **kwds):
        """
        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]

        We create a skew polynomial from a list::

            sage: S([t,1])
            x + t

        from another skew polynomial::

            sage: S(x^2 + t)
            x^2 + t

        from a constant::

            sage: x = S(t^2 + 1); x
            t^2 + 1
            sage: x.parent() is S
            True
        """
        SkewPolynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self.__coeffs = []
            return

        R = parent.base_ring()
#        if PY_TYPE_CHECK(x, list):
        if type(x) is list:
            if check:
                self.__coeffs = [R(t) for t in x]
                self.__normalize()
            else:
                self.__coeffs = x
            return

#        if PY_TYPE_CHECK(x, SkewPolynomial):
        if type(x) is SkewPolynomial:
            if (<Element>x)._parent is self._parent:
                x = list(x.list())
            elif R.has_coerce_map_from((<Element>x)._parent):# is R or (<Element>x)._parent == R:
                try:
                    if x.is_zero():
                        self.__coeffs = []
                        return
                except (AttributeError, TypeError):
                    pass
                x = [x]
            else:
                self.__coeffs = [R(a, **kwds) for a in x.list()]
                if check:
                    self.__normalize()
                return

#        elif PY_TYPE_CHECK(x, int) and x == 0:
        elif type(x) is int and x == 0:
            self.__coeffs = []
            return

        elif isinstance(x, dict):
            x = self._dict_to_list(x, R.zero())
#            x = self._dict_to_list(x, R.zero_element())

        elif not isinstance(x, list):
            # We trust that the element constructors do not send x=0
#            if x:
            x = [x] # constant polynomials
#            else:
#                x = [] # zero polynomial
        if check:
            self.__coeffs = [R(z, **kwds) for z in x]
            self.__normalize()
        else:
            self.__coeffs = x


    cdef list _list_c(self):
        """
        Return the list of the underlying elements of self.

        .. WARNING::

            It is a priori not a copy; do not modify this list!
        """
        return self.__coeffs


    cdef SkewPolynomial _new_c(self, list coeffs, Parent P, char check=0):
        """
        Fast creation of a new skew polynomial
        """
#        cdef SkewPolynomial_generic_dense f = <SkewPolynomial_generic_dense>PY_NEW_SAME_TYPE(self)
        cdef type t = type(self)
        cdef SkewPolynomial_generic_dense f = t.__new__(t)
        f._parent = P
        f.__coeffs = coeffs
        if check:
            f.__normalize()
        return f


    cdef void __normalize(self):
        x = self.__coeffs
        cdef Py_ssize_t n = len(x) - 1
        while n >= 0 and not x[n]:
            del x[n]
            n -= 1


    # Basic operations in place
    # -------------------------

    cdef void _inplace_add(self, SkewPolynomial_generic_dense right):
        """
        Replace self by self+right (only for internal use).
        """
        cdef Py_ssize_t i, min
        x = (<SkewPolynomial_generic_dense>self).__coeffs
        y = (<SkewPolynomial_generic_dense>right).__coeffs
        if len(x) > len(y):
            for i from 0 <= i < len(y):
                x[i] += y[i]
        else:
            x += y[len(x):]
            for i from 0 <= i < len(x):
                x[i] += y[i]
            x += y[len(x):]
        if len(x) == len(y):
            self.__normalize()

    cdef void _inplace_sub(self, SkewPolynomial_generic_dense right):
        """
        Replace self by self-right (only for internal use).
        """
        cdef Py_ssize_t i, min
        cdef list x = (<SkewPolynomial_generic_dense>self).__coeffs
        cdef list y = (<SkewPolynomial_generic_dense>right).__coeffs
        if len(x) >= len(y):
            for i from 0 <= i < len(y):
                x[i] -= y[i]
        else:
            for i from 0 <= i < len(x):
                x[i] -= y[i]
            x += [-c for c in y[len(x):]]
        if len(x) == len(y):
            self.__normalize()

    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right):
        """
        Replace self by self*right (only for internal use).
        """
        cdef list x = (<SkewPolynomial_generic_dense>self).__coeffs
        cdef list y = (<SkewPolynomial_generic_dense>right).__coeffs
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t d1 = len(x)-1, d2 = len(y)-1
        parent = self._parent
        if d2 == -1:
            (<SkewPolynomial_generic_dense>self).__coeffs = [ ]
        elif d1 >= 0:
            for k from d1 < k <= d1+d2:
                start = 0 if k <= d2 else k-d2 # max(0, k-d2)
                sum = x[start] * parent.twist_map(start)(y[k-start])
                for i from start < i <= d1:
                    sum += x[i] * parent.twist_map(i)(y[k-i])
                x.append(sum)
            for k from d1 >= k >= 0:
                start = 0 if k <= d2 else k-d2 # max(0, k-d2)
                end = k if k <= d1 else d1 # min(k, d1)
                sum = x[start] * parent.twist_map(start)(y[k-start])
                for i from start < i <= end:
                    sum += x[i] * parent.twist_map(i)(y[k-i])
                x[k] = sum

    cdef void _inplace_lmul(self, SkewPolynomial_generic_dense left):
        """
        Replace self by left*self (only for internal use).
        """
        cdef list x = (<SkewPolynomial_generic_dense>self).__coeffs
        cdef list y = (<SkewPolynomial_generic_dense>left).__coeffs
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t d1 = len(x)-1, d2 = len(y)-1
        parent = self._parent
        if d2 == -1:
            (<SkewPolynomial_generic_dense>self).__coeffs = [ ]
        elif d1 >= 0:
            for k from d1 < k <= d1+d2:
                start = 0 if k <= d2 else k-d2 # max(0, k-d2)
                sum = parent.twist_map(k-start)(x[start]) * y[k-start]
                for i from start < i <= d1:
                    sum += parent.twist_map(k-i)(x[i]) * y[k-i]
                x.append(sum)
            for k from d1 >= k >= 0:
                start = 0 if k <= d2 else k-d2 # max(0, k-d2)
                end = k if k <= d1 else d1 # min(k, d1)
                sum = parent.twist_map(k-start)(x[start]) * y[k-start]
                for i from start < i <= end:
                    sum += parent.twist_map(k-i)(x[i]) * y[k-i]
                x[k] = sum


    # Fast exponentiation
    # -------------------

    cdef void _inplace_pow(self, Py_ssize_t n):
        """
        Replace self by self**n.
        """
        while n & 1 == 0:
            self._inplace_rmul(self)
            n = n >> 1
        cdef SkewPolynomial_generic_dense selfpow = <SkewPolynomial_generic_dense>self._new_c(list(self.__coeffs),self._parent)
        n = n >> 1
        while n != 0:
            selfpow._inplace_rmul(selfpow)
            if n&1 == 1:
                self._inplace_rmul(selfpow)
            n = n >> 1


    cpdef _pow_(self,exp,modulus=None,side=Right):
        """
        INPUT:

        -  ``exp`` -- an Integer

        -  ``modulus`` -- a skew polynomial over the same ring (default: None)

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        If ``modulus`` is None, return ``self**exp``.

        Otherwise, return the remainder of self**exp in the ``side``
        euclidean division by ``modulus``.

        REMARK:

        The quotient of the underlying skew polynomial ring by the
        principal ideal generated by ``modulus`` is in general *not*
        a ring.

        As a consequence, Sage first computes exactly ``self**exp``
        and then reduce it modulo ``modulus``.

        However, if the base ring is a finite field, Sage uses the
        following optimized algorithm:

        #. One first compute a central skew polynomial `N` which is
           divisible by ``modulus``. (Since `N` lies in center, the
           quotient `K[X,\sigma]/N` inherits a ring structure.)

        #. One compute ``self**exp`` in the quotient ring `K[X,\sigma]/N`

        #. One reduce modulo ``modulus`` the result computed in the
           previous step

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: b = a^10  # short form for ``a._pow_(10)``
            sage: b == a*a*a*a*a*a*a*a*a*a
            True

            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: br = a._pow_(10,modulus); br
            (t^2 + t)*x^2 + (3*t^2 + 1)*x + t^2 + t
            sage: br == b.rem(modulus)
            True
            sage: bl = a._pow_(10,modulus,side=Left); bl
            (4*t^2 + 2*t + 3)*x^2 + (3*t^2 + 1)*x + 2*t + 3
            sage: bl == b.rem(modulus,side=Left)
            True

            sage: a._pow_(10^100,modulus)  # quite fast
            (3*t^2 + 3)*x^2 + (t^2 + 2*t + 4)*x + 4*t^2 + 2*t + 1
        """
        cdef SkewPolynomial_generic_dense r

#        if not PY_TYPE_CHECK_EXACT(exp, Integer) or \
#                PY_TYPE_CHECK_EXACT(exp, int):
        if not type(exp) is Integer or \
                type(exp) is int:
                    try:
                        exp = Integer(exp)
                    except TypeError:
                        raise TypeError("non-integral exponents not supported")

        if self.degree() <= 0:
            return self.parent()(self[0]**exp)
        if exp == 0:
            return self.parent()(1)
        if exp < 0:
            return (~self).pow(-exp,modulus,side=side)

        if self == self.parent().gen(): # special case x**n should be faster!
            P = self.parent()
            R = P.base_ring()
#            v = [R.zero_element()]*exp + [R.one_element()]
            v = [R.zero()]*exp + [R.one()]
            r = <SkewPolynomial_generic_dense>self._parent(v)
        else:
            r = <SkewPolynomial_generic_dense>self._new_c(list(self.__coeffs),self._parent)
            sig_on()
            r._inplace_pow(exp)
            sig_off()

        if modulus:
            sig_on()
            if side == Right:
                r._inplace_rrem(modulus)
            else:
                r._inplace_lrem(modulus)
            sig_off()

        return r


    def __pow__(self,exp,modulus):
        """
        INPUT:

        -  ``exp`` -- an Integer

        -  ``modulus`` -- a skew polynomial over the same ring
           (default: None)

        OUTPUT:

        If ``modulus`` is None, return ``self**exp``.

        Otherwise, return the remainder of self**exp in the right
        euclidean division by ``modulus``.

        .. SEEALSO::

            :meth:`~sage.rings.polynomial.skew_polynomial_element._pow_`

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: b = a^10
            sage: b == a*a*a*a*a*a*a*a*a*a
            True

            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: bmod = a._pow_(10,modulus); bmod
            (t^2 + t)*x^2 + (3*t^2 + 1)*x + t^2 + t
            sage: bmod == b.rem(modulus)
            True
        """
        return self._pow_(exp,modulus)


def make_generic_skew_polynomial(parent, coeffs):
    return parent(coeffs)


cdef class ConstantSkewPolynomialSection(Map):
    cpdef Element _call_(self, x):
        if x.degree() <= 0:
            try:
                return <Element>(x.constant_coefficient())
            except AttributeError:
                return <Element>((<SkewPolynomial>x).constant_coefficient())
        else:
            raise TypeError("not a constant polynomial")


cdef class SkewPolynomialBaseringInjection(RingHomomorphism):
    def __init__(self, domain, codomain):
        assert codomain.base_ring() is domain, "domain must be basering"
        RingHomomorphism.__init__(self, Hom(domain,codomain))
        self._an_element = codomain.gen()
        self._repr_type_str = "Polynomial base injection"
        self._new_constant_poly_ = self._an_element._new_constant_poly

    def an_element(self):
        return self._an_element

    def new_constant_poly_(self):
        return self._new_constant_poly_

    cpdef Element _call_(self, x):
        return self._new_constant_poly_(x, self._codomain)
#        return self._new_constant_poly_(x, self._codomain)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            return self._codomain._element_constructor_(x, *args, **kwds)
        except AttributeError:
            # if there is no element constructor, there is a custom call method.
            return self._codomain(x, *args, **kwds)

    def section(self):
        return ConstantSkewPolynomialSection(self._codomain, self._domain)
