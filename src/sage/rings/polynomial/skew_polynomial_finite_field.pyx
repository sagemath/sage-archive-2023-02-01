r"""
This module implements skew polynomials over finite fields.

Let `k` be a finite field and `\sigma` be a ring automorphism
of `k` (i.e. a power of the Frobenius endomorphism). Let
Put `S = k[X,\sigma]`: as an addtive group, it is the usual ring
of polynomials with coefficients in `k` and the multiplication
on `S` is defined by the rule `X * a = \sigma(a) * X`.

We recall that:

#. `S` is a left (resp. right) euclidean noncommutative ring

#. in particular, every left (resp. right) ideal is principal

Since `k` is a finite field, we have the additional following
properties:

#. the center of `S`, denoted by `Z`, is the univariate
   polynomial ring over `k^\sigma` (`k` fixed by `\sigma`)
   in the variable `X^r` where `r` is the order of `\sigma`

#. `S[1/X]` is an Azumaya algebra over `Z[1/X^r]` (i.e. etale
   locally `S[1/X]` is a matrix algebra over `Z[1/X^r]`)

#. in particular, we have a reduced norm map `N` from `S[1/X]`
   to `Z[1/X^r]` (etale locally, it is the determinant); one
   can prove that it maps `S` to `Z`

#. `N` has very good properties regarding to factorizations; in
   particular:

   #. if `a` is a skew polynomial, `a` always divides `N(a)`

   #. if `a` is a skew polynomial, any factorization of `N(a)`
      (in any order) lifts to a factorization of `a` (and we
      have a precise control on the number of such lifts); as
      a consequence there is an explicit (but complicated)
      formula counting the number of factorizations of a skew
      polynomial.


.. TODO::

    Try to replace as possible ``finite field`` by ``field
    endowed with a finite order twist morphism``. It may cause
    new phenomena due to the non trivality of the Brauer group.

EXAMPLES::

We illustrate some properties listed above::

    sage: k.<t> = GF(5^3)
    sage: Frob = k.frobenius_endomorphism()
    sage: S.<x> = k['x',Frob]; S
    Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
    sage: Z = S.center(); Z
    Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
    Univariate Polynomial Ring in (x^3) over Finite Field of size 5
    sage: a = x^5 + (2*t^2 + t + 1)*x^4 + (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2

We compute the reduced norm of `a`::

    sage: N = a.reduced_norm(); N
    (x^3)^5 + 3*(x^3)^4 + 2*(x^3)^3 + 3*(x^3) + 4

Note that the parent of `N` is the center `Z` (and not `S`)::

    sage: N.parent()
    Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
    Univariate Polynomial Ring in (x^3) over Finite Field of size 5
    sage: N.parent() == Z
    True

    sage: S(N)  # coercion of N into S
    x^15 + 3*x^12 + 2*x^9 + 3*x^3 + 4

`N` is a mutiple of `a`::

    sage: S(N).is_divisible_by(a)
    True
    sage: S(N).is_divisible_by(a,side=Left)
    True

.. NOTE::

    We really need to coerce first `N` into `S`. Otherwise an
    error occurs::

        sage: N.is_divisible_by(a)
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.polynomial.skew_polynomial_element.CenterSkewPolynomial_generic_dense' object has no attribute 'is_divisible_by'

As a polynomial in `x^r` (here `r = 3`), ``N`` factors as a produit
of two irreducible polynomials::

    sage: N.factor()
    ((x^3)^2 + (x^3) + 1) * ((x^3)^3 + 2*(x^3)^2 + 4*(x^3) + 4)

And so does `a`::

    sage: F = a.factor(); F
    (x^3 + (t^2 + 2*t + 1)*x^2 + (4*t^2 + t + 2)*x + 3*t^2 + t + 4) * (x^2 + (t^2 + 4*t)*x + 2*t)

We can check that each of these two factors of `a` corresponds
to a factor of `N`::

    sage: F[0][0].reduced_norm()
    (x^3)^3 + 2*(x^3)^2 + 4*(x^3) + 4
    sage: F[1][0].reduced_norm()
    (x^3)^2 + (x^3) + 1

Actually, `a` has exactly two factorizations corresponding to the
two possible orderings of the irreducible factors of `N`::

    sage: a.count_factorizations()
    2
    sage: for F in a.factorizations(): print F
    (x^3 + (t^2 + 2*t + 1)*x^2 + (4*t^2 + t + 2)*x + 3*t^2 + t + 4) * (x^2 + (t^2 + 4*t)*x + 2*t)
    (x^2 + (2*t + 2)*x + t + 1) * (x^3 + (3*t^2 + 1)*x^2 + (2*t + 4)*x + 4*t^2 + 3*t + 4)

AUTHOR::

- Xavier Caruso (2012-06-29)
"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************


include "../../ext/stdsage.pxi"
#include "../../ext/interrupt.pxi"
include "cysignals/signals.pxi"

import operator, copy, re

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.functions.other import ceil

from sage.structure.parent cimport Parent

from skew_polynomial_element cimport SkewPolynomial_generic_dense
from sage.matrix.constructor import matrix
from sage.matrix.constructor import zero_matrix
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.factorization import Factorization

from skew_polynomial_element cimport SkewPolynomial
from polynomial_ring import PolynomialRing_general
from polynomial_ring_constructor import PolynomialRing

from sage.rings.ring cimport Ring
from sage.structure.element cimport RingElement

from sage.structure.side import Left, Right
from sage.rings.infinity import Infinity

from sage.combinat.q_analogues import q_jordan
from sage.functions.other import factorial
from sage.combinat.permutation import Permutations
from sage.combinat.partition import Partition
from sage.misc.mrange import xmrange_iter


cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_generic_dense):

    def __init__(self, parent, x=None, int check=1, is_gen=False, int construct=0, **kwds):
        SkewPolynomial_generic_dense.__init__ (self, parent, x, check, is_gen, construct, **kwds)
        self._init_cache()


    cdef inline void _init_cache(self):
        """
        Initialize cached variables (set them to None).
        """
        self._conjugates = [ self.__coeffs ]
        self._norm = None
        self._norm_factor = None
        self._optbound = None
        self._rdivisors = None
        self._types = None
        self._factorization = None


    cdef SkewPolynomial _new_c(self, list coeffs, Parent P, char check=0):
        """
        Fast creation of a new skew polynomial
        """
#        cdef SkewPolynomial_finite_field_dense f = <SkewPolynomial_finite_field_dense>PY_NEW_SAME_TYPE(self)
        cdef type t = type(self)
        cdef SkewPolynomial_finite_field_dense f = t.__new__(t)
        f._parent = P
        f.__coeffs = coeffs
        f._init_cache()
        if check:
            f.__normalize()
        return f
    
    def norm(self):
        """
        Return norm of `self`
        """
        return self._norm

    def norm_factor(self):
        """
        Return norm_factor of `self`
        """
        return self._norm_factor

    # Skew multiplication
    # -------------------

    def _mul_karatsuba(self,right,cutoff=None):
        """
        Karatsuba multiplication

        INPUT:

        - ``right`` -- an other skew polynomial in the same ring

        - ``cutoff`` -- ``None``, an integer or Infinity (default: None)

        .. WARNING::

            ``cutoff`` need to be greater than or equal to the order of the
            twist map acting on the base ring of the underlying skew polynomial
            ring.

        OUTPUT:

        The result of the product self*right (computed by a variant of
        Karatsuba`s algorithm)

        .. NOTE::

            if ``cutoff`` is None, use the default cutoff which is the
            maximum between 150 and the order of the twist map.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = S.random_element(5000)
            sage: b = S.random_element(5000)
            sage: timeit("c = a._mul_karatsuba(b)")  # random, long time
            5 loops, best of 3: 659 ms per loop
            sage: timeit("c = a._mul_classical(b)")  # random, long time
            5 loops, best of 3: 1.9 s per loop
            sage: a._mul_karatsuba(b) == a._mul_classical(b)
            True

        The operator ``*`` performs Karatsuba multiplication::

            sage: timeit("c = a*b")  # random, long time
            5 loops, best of 3: 653 ms per loop
        """
        karatsuba_class = self._parent._karatsuba_class
        if cutoff != None:
            save_cutoff = karatsuba_class.get_cutoff()
            karatsuba_class.set_cutoff(cutoff)
        res = karatsuba_class.mul(self,right)
        if cutoff != None:
            karatsuba_class.set_cutoff(save_cutoff)
        return res


    def _mul_karatsuba_matrix(self,right):
        """
        Karatsuba multiplication with multiplication step based
        on an isomorphism between a quotient of the underlying
        skew polynomial ring and a ring of matrices.

        INPUT:

        - ``right`` -- an other skew polynomial in the same ring

        OUTPUT:

        The result of the product self*right

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = S.random_element(degree=20)
            sage: b = S.random_element(degree=20)
            sage: a._mul_karatsuba_matrix(b) == a*b
            True

        This routine is only efficient when the twisting map (here
        ``Frob``) has a large order `r` and the degrees of ``self``
        and ``other`` have a very special shape (just below a power
        of `2` times `r * floor(r/2)`)::

            sage: k.<t> = GF(5^40)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = S.random_element(degree=799)
            sage: b = S.random_element(degree=799)
            sage: timeit("c = a*b")  # random, long time
            5 loops, best of 3: 23.1 s per loop
            sage: timeit("c = a._mul_karatsuba_matrix(b)")  # random, long time
            5 loops, best of 3: 12.2 s per loop
        """
        karatsuba_class = self._parent._karatsuba_class
        return karatsuba_class.mul_matrix(self,right)


    cpdef SkewPolynomial_finite_field_dense _mul_central(self, SkewPolynomial_finite_field_dense right):
        r"""
        Return self * right

        .. WARNING::

            Do you use this function! It is very slow due to a quite
            slow interface with ``polynomial_zz_pex``.

        ALGORITHM::

        Notations::

        -  `S` is the underlyling skew polynomial ring

        -  `x` is the variable on `S`

        -  `k` is the base ring of `S` (it is a finite field)

        -  `\sigma` is the twisting automorphism acting on `k`

        -  `r` is the order of `\sigma`

        -  `t` is a generator of `k` over `k^\sigma`

        #. We decompose the polynomial ``right`` as follows::

           .. MATH::

               right = \sum_{i=0}^{r-1} \sum_{j=0}^{r-1} y_{i,j} t^j x^i

           where `y_{i,j}` are polynomials in the center `k^\sigma[x^r]`.

        #. We compute all products `z_{i,j} = left * y_{i,j}`; since
           all `y_{i,j}` lie in the center, we can compute all these
           products as if `left` was a commutative polynomial (and we
           can therefore use fast algorithms like FFT and/or fast
           implementations)

        #. We compute and return the sum

           .. MATH::

               \sum_{i=0}^{r-1} \sum_{j=0}^{r-1} z_{i,j} t^j x^i

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = S.random_element(degree=10)
            sage: b = S.random_element(degree=10)
            sage: a._mul_central(b) == a*b
            True

        TESTS::

        Here is an example where `k^\sigma` is not a prime field::

            sage: k.<t> = GF(5^6)
            sage: Frob = k.frobenius_endomorphism(2)
            sage: S.<x> = k['x',Frob]
            sage: a = S.random_element(degree=10)
            sage: b = S.random_element(degree=10)
            sage: a._mul_central(b) == a*b
            True
        """
        skew_ring = self._parent
        base_ring = skew_ring.base_ring()
        commutative_ring = PolynomialRing(skew_ring.base_ring(),name='x')
        cdef RingElement c
        cdef RingElement zero = base_ring(0)
        cdef Py_ssize_t i, j, k
        cdef Py_ssize_t order = skew_ring._order
        cdef Py_ssize_t degree = base_ring.degree()

        left = commutative_ring(self.__coeffs)
        cdef list y = [ c.polynomial() for c in right.__coeffs ]
        cdef Py_ssize_t leny = len(y)
        cdef list yc = leny * [zero]
        cdef list res = (leny + len(self.__coeffs) - 1) * [zero]
        cdef list term
        cdef list twist = [ base_ring.gen() ]
        for i from 0 <= i < order-1:
            twist.append(skew_ring.twist_map(1)(twist[i]))
        for i from 0 <= i < order:
            for j from 0 <= j < degree:
                for k from i <= k < leny by order:
                    yc[k] = y[k][j]
                term = (left * commutative_ring(yc)).list()
                for k from i <= k < len(term):
                    res[k] += term[k] * twist[(k-i)%order]**j
            for k from i <= k < leny by order:
                yc[k] = zero
        return self._new_c(res,skew_ring,1)


    cpdef RingElement _mul_(self, RingElement right):
        """
        Compute self * right (in this order)

        .. NOTE::

            Use skew Karatsuba's algorithm for skew
            polynomials of large degrees.

        INPUT:

        -  right -- a skew polynomial in the same ring

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t; a
            x^2 + t
            sage: b = x^2 + (t + 1)*x; b
            x^2 + (t + 1)*x
            sage: a * b
            x^4 + (3*t^2 + 2)*x^3 + t*x^2 + (t^2 + t)*x
            sage: a * b == b * a
            False
        """
        return self._parent._karatsuba_class.mul(self,right)


    def _mul_classical(self,right):
        """
        Compute self * right (in this order) using the
        skew SchoolBook algorithm.

        INPUT:

        -  ``right`` -- a skew polynomial in the same ring

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t; a
            x^2 + t
            sage: b = x^2 + (t + 1)*x; b
            x^2 + (t + 1)*x
            sage: a._mul_classical(b)
            x^4 + (3*t^2 + 2)*x^3 + t*x^2 + (t^2 + t)*x
            sage: a * b == b * a
            False
        """
        karatsuba_class = self._parent._karatsuba_class
        save_cutoff = karatsuba_class.get_cutoff()
        karatsuba_class.set_cutoff(Infinity)
        res = karatsuba_class.mul(self,right)
        karatsuba_class.set_cutoff(save_cutoff)
        return res


    cpdef rquo_rem_karatsuba(self, RingElement other, cutoff=None):
        """
        Right euclidean division based on Karatsuba's algorithm.

        DO NOT USE THIS! It is not efficient for usual degrees!

        .. TODO::

            Try to understand why...

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: a = S.random_element(2000)
            sage: b = S.random_element(1000)
            sage: timeit("q,r = a.rquo_rem_karatsuba(b)")  # random, long time
            5 loops, best of 3: 104 ms per loop
            sage: timeit("q,r = a.rquo_rem(b)")  # random, long time
            5 loops, best of 3: 79.6 ms per loop
            sage: a.rquo_rem(b) == a.rquo_rem_karatsuba(b)
            True

            sage: a = S.random_element(10000)
            sage: b = S.random_element(5000)
            sage: timeit("q,r = a.rquo_rem_karatsuba(b)")  # random, long time
            5 loops, best of 3: 1.79 s per loop
            sage: timeit("q,r = a.rquo_rem(b)")  # random, long time
            5 loops, best of 3: 1.93 s per loop
            sage: a.rquo_rem(b) == a.rquo_rem_karatsuba(b)
            True
        """
        karatsuba_class = self._parent._karatsuba_class
        if cutoff != None:
            save_cutoff = karatsuba_class.get_cutoff()
            karatsuba_class.set_cutoff(cutoff)
        res = karatsuba_class.div(self,other)
        if cutoff != None:
            karatsuba_class.set_cutoff(save_cutoff)
        return res


    # We improve some functions
    # -------------------------

    def rquo_rem(self,other):
        """
        DEFINITION:

        Let `a` and `b` be two skew polynomials over the same
        ring. The *right euclidean division* of `a` by `b` is
        a couple `(q,r)` such that

        -  `a = q*b + r`

        -  the degree of `r` is less than the degree of `b`

        `q` (resp. `r`) is called the *quotient* (resp. the
        remainder) of this euclidean division.

        If the leading coefficient of `b` is a unit (e.g. if
        `b` is monic) then `q` and `r` exist and are unique.

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
        cdef list a = self.list()
        cdef list b = other.list()
        cdef Py_ssize_t i, j
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            return parent(0), self
        cdef RingElement inv = ~b[db]
        cdef list q = [ ]
        cdef Py_ssize_t order = parent._order
        cdef list twinv = [ inv ], twb = [ b ]
        cdef RingElement c, x
        for i from 0 <= i < min(da-db,order-1):
            twinv.append(parent.twist_map()(twinv[i]))
            twb.append([ parent.twist_map()(x) for x in twb[i] ])
        for i from da-db >= i >= 0:
            c = twinv[i%order] * a[i+db]
            for j from 0 <= j < db:
                a[i+j] -= c * twb[i%order][j]
            q.append(c)
        q.reverse()
        return parent(q), parent(a[:db])


    cdef SkewPolynomial_finite_field_dense _rgcd(self,SkewPolynomial_finite_field_dense other):
        """
        Fast right gcd.
        """
        cdef SkewPolynomial_finite_field_dense A = self
        cdef SkewPolynomial_finite_field_dense B = other
        cdef SkewPolynomial_finite_field_dense swap
        if len(B.__coeffs):
            A = <SkewPolynomial_finite_field_dense>self._new_c(A.__coeffs[:],A._parent)
            B = <SkewPolynomial_finite_field_dense>B._new_c(B.__coeffs[:],B._parent)
            while len(B.__coeffs):
                A._inplace_rrem(B)
                swap = A; A = B; B = swap
            return A
        else:
            return self


    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right):
        """
        Replace self by self*right
        """
        self.__coeffs = self._parent._karatsuba_class.mul_list(self.__coeffs,right.__coeffs)
        self._init_cache()


    cdef void _inplace_lmul(self, SkewPolynomial_generic_dense left):
        """
        Replace self by left*self
        """
        self.__coeffs = self._parent._karatsuba_class.mul_list(left.__coeffs,self.__coeffs)
        self._init_cache()


    cpdef _pow_(self,right,modulus=None,side=Right):
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

            sage: a._pow_(10^100,modulus)  # rather fast
            (3*t^2 + 3)*x^2 + (t^2 + 2*t + 4)*x + 4*t^2 + 2*t + 1
        """
        cdef SkewPolynomial_finite_field_dense r

        if not isinstance(right, Integer) or isinstance(right, int):
            try:
                right = Integer(right)
            except TypeError:
                raise TypeError("non-integral exponents not supported")

        if self.degree() <= 0:
            return self.parent()(self[0]**right)
        if right == 0:
            return self.parent()(1)
        if right < 0:
            return (~self).pow(-right,modulus,side=side)

        if self == self.parent().gen(): # special case x**n should be faster!
            P = self.parent()
            R = P.base_ring()
            v = [R.zero()]*right + [R.one()]
#            v = [R.zero_element()]*right + [R.one_element()]
            r = <SkewPolynomial_generic_dense>self._parent(v)
            sig_on()
            if modulus:
                if side is Right:
                    r._inplace_rrem(modulus)
                else:
                    r._inplace_lrem(modulus)
                r._init_cache()
            sig_off()
            return r

        mod = modulus
        if not modulus is None:
            try:
                mod = self.parent()(mod.bound())
            except NotImplementedError:
                mod = None
        r = <SkewPolynomial_generic_dense>self._new_c(self.__coeffs,self._parent)
        sig_on()
        if mod:
            r._inplace_pow_mod(right,mod)
        else:
            r._inplace_pow(right)
        if (not modulus is None) and modulus != mod:
            if side is Right:
                r._inplace_rrem(modulus)
            else:
                r._inplace_lrem(modulus)
            r._init_cache()
        sig_off()
        return r


    # Inplace functions
    # -----------------

    #cdef void _inplace_conjugate(self,n):
    #    cdef Py_ssize_t i
    #    cdef Morphism twist = <Morphism>self._parent.twist_map(n)
    #    cdef RingElement x
    #    for i from 0 <= i < len(self.__coeffs):
    #        x = twist(self.__coeffs[i])
    #        self.__coeffs[i] = x


    cdef void _inplace_lrem(self, SkewPolynomial_finite_field_dense other):
        """
        Replace self by the remainder in the left euclidean division
        of self by other (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other).__coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j
        cdef RingElement c, inv
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da >= db:
            inv = ~b[db]
            for i from da-db >= i >= 0:
                c = parent.twist_map(-db)(inv*a[i+db])
                for j from 0 <= j < db:
                    a[i+j] -= b[j] * parent.twist_map(j)(c)
            del a[db:]
            self.__normalize()
        self._init_cache()


    cdef void _inplace_rrem(self, SkewPolynomial_finite_field_dense other):
        """
        Replace self by the remainder in the right euclidean division
        of self by other (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other).__coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j, order
        cdef RingElement c, x, inv
        cdef list twinv, twb
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da >= db:
            order = parent._order
            inv = ~b[db]
            twinv = [ inv ]
            for i from 0 <= i < min(da-db,order-1):
                twinv.append(parent.twist_map()(twinv[i]))
            twb = (<SkewPolynomial_finite_field_dense>other)._conjugates
            for i from len(twb)-1 <= i < min(da-db,order-1):
                twb.append([ parent.twist_map()(x) for x in twb[i] ])
            for i from da-db >= i >= 0:
                c = twinv[i%order] * a[i+db]
                for j from 0 <= j < db:
                    a[i+j] -= c * twb[i%order][j]
            del a[db:]
            self.__normalize()
        self._init_cache()


    cdef void _inplace_lfloordiv(self, SkewPolynomial_finite_field_dense other):
        """
        Replace self by the quotient in the left euclidean division
        of self by other (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other).__coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j, deb
        cdef RingElement c, inv
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            (<SkewPolynomial_finite_field_dense>self).__coeffs = [ ]
        else:
            inv = ~b[db]
            for i from da-db >= i >= 0:
                c = a[i+db] = parent.twist_map(-db)(inv*a[i+db])
                if i < db: deb = db
                else: deb = i
                for j from deb <= j < db+i:
                    a[j] -= b[j-i] * parent.twist_map(j-i)(c)
            del a[:db]
            self.__normalize()
        self._init_cache()


    cdef void _inplace_rfloordiv(self, SkewPolynomial_finite_field_dense other):
        """
        Replace self by the quotient in the right euclidean division
        of self by other (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other).__coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j, deb, order
        cdef RingElement c, x, inv
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            (<SkewPolynomial_finite_field_dense>self).__coeffs = [ ]
        else:
            order = parent._order
            inv = ~b[db]
            twinv = [ inv ]
            for i from 0 <= i < min(da-db,order-1):
                twinv.append(parent.twist_map()(twinv[i]))
            twb = (<SkewPolynomial_finite_field_dense>other)._conjugates
            for i from len(twb)-1 <= i < min(da-db,order-1):
                twb.append([ parent.twist_map()(x) for x in twb[i] ])
            for i from da-db >= i >= 0:
                c = a[i+db] = twinv[i%order] * a[i+db]
                if i < db: deb = db
                else: deb = i
                for j from deb <= j < db+i:
                    a[j] -= c * twb[i%order][j-i]
            del a[:db]
            self.__normalize()
        self._init_cache()


    cdef void _inplace_pow_mod(self, Integer n, SkewPolynomial_finite_field_dense mod):
        """
        Replace self by the remainder in the euclidean
        division of self**n by mod (only for internal use).

        .. WARNING::

            Assume that mod is central
        """
        while n&1 == 0:
            self._inplace_rmul(self)
            self._inplace_rrem(mod)
            n = n >> 1
        cdef SkewPolynomial_finite_field_dense selfpow = <SkewPolynomial_finite_field_dense>self._new_c(self.__coeffs[:], self._parent)
        n = n >> 1
        while n != 0:
            selfpow._inplace_rmul(selfpow)
            selfpow._inplace_rrem(mod)
            if n&1 == 1:
                self._inplace_rmul(selfpow)
                self._inplace_rrem(mod)
            n = n >> 1


    cdef void _inplace_lmonic(self):
        """
        Replace self by ``self.lmonic()`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef Py_ssize_t da = len(a)-1, i
        cdef RingElement inv = ~a[da]
        parent = self._parent
        a[da] = parent.base_ring()(1)
        for i from 0 <= i < da:
            a[i] *= parent.twist_map(i-da)(inv)
        self._init_cache()


    cdef void _inplace_rmonic(self):
        """
        Replace self by ``self.rmonic()`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef Py_ssize_t da = len(a)-1, i
        cdef RingElement inv = ~a[da]
        a[da] = self._parent.base_ring()(1)
        for i from 0 <= i < da:
            a[i] *= inv
        self._init_cache()


    cdef void _inplace_rgcd(self,SkewPolynomial_finite_field_dense other):
        """
        Replace self by its right gcd with other (only for internal use).
        """
        cdef SkewPolynomial_finite_field_dense B
        cdef list swap
        if len(other.__coeffs):
            B = <SkewPolynomial_finite_field_dense>self._new_c(other.__coeffs[:],other._parent)
            while len(B.__coeffs):
                B._conjugates = [ B.__coeffs ]
                self._inplace_rrem(B)
                swap = self.__coeffs
                self.__coeffs = B.__coeffs
                B.__coeffs = swap
        self._init_cache()


    cdef SkewPolynomial_finite_field_dense _rquo_inplace_rem(self, SkewPolynomial_finite_field_dense other):
        """
        Replace self by the remainder in the right euclidean division
        of self by other and return the quotient (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other).__coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j
        cdef RingElement c, inv
        cdef list q
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            return self._new_c([],self._parent)
        inv = ~b[db]
        q = [ ]
        for i from da-db >= i >= 0:
            c = parent.twist_map(i)(inv) * a[i+db]
            q.append(c)
            for j from 0 <= j < db:
                a[i+j] -= c * parent.twist_map(i)(b[j])
        del a[db:]
        self.__normalize()
        self._init_cache()
        q.reverse()
        return self._new_c(q,self._parent)


    cdef Py_ssize_t _val_inplace_unit(self):
        """
        Return `v` the valuation of self and replace self by
        self >> v (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self).__coeffs
        cdef Py_ssize_t val = 0
        if len(a) < 0:
            return -1
        while a[0].is_zero():
            del a[0]
            val += 1
        self._init_cache()
        return val


    # lowest common divisor: other functions
    # --------------------------------------

    cdef SkewPolynomial_finite_field_dense _coeff_llcm(self, SkewPolynomial_finite_field_dense other):
        """
        (only for internal use)
        """
        cdef SkewPolynomial_finite_field_dense A = <SkewPolynomial_finite_field_dense>self
        cdef SkewPolynomial_finite_field_dense B = <SkewPolynomial_finite_field_dense>other
        if not len(A.__coeffs) and not len(B.__coeffs):
            raise ZeroDivisionError
        A = <SkewPolynomial_finite_field_dense>A._new_c(A.__coeffs[:],A._parent)
        B = <SkewPolynomial_finite_field_dense>A._new_c(B.__coeffs[:],B._parent)
        cdef parent = self._parent
        cdef RingElement one = self.base_ring().one()
#        cdef RingElement one = self.base_ring().one_element()
        cdef SkewPolynomial_finite_field_dense U = self._new_c([one],parent)
        cdef SkewPolynomial_finite_field_dense V = self._new_c([],parent)
        cdef SkewPolynomial_finite_field_dense Q, R
        cdef list swap
        while len(B.__coeffs):
            Q = A._rquo_inplace_rem(B)
            Q._inplace_rmul(V)
            U._inplace_sub(Q)
            swap = A.__coeffs
            A.__coeffs = B.__coeffs
            B.__coeffs = swap
            swap = U.__coeffs
            U.__coeffs = V.__coeffs
            V.__coeffs = swap
        return V


    # Specific functions
    # ------------------

    cdef Matrix_dense _matphir_c(self):
        r"""
        Return the matrix of the multiplication by `X^r` on
        the quotient `K[X,\sigma] / K[X,\sigma]*self`.
        """
        cdef Py_ssize_t i, j, col, exp, n
        cdef Py_ssize_t d = self.degree()
        cdef Py_ssize_t r = self.parent()._order
        parent = self._parent
        cdef k = parent.base_ring()
        cdef Matrix_dense phir = <Matrix_dense?>zero_matrix(k,d,d)
        cdef RingElement zero = k(0)
        cdef RingElement one = k(1)
        if r < d:
            for i from 0 <= i < d-r:
                phir.set_unsafe(r+i,i,one)
            col = d-r
            exp = d
        else:
            col = 0
            exp = r
        cdef SkewPolynomial_finite_field_dense powx = <SkewPolynomial_finite_field_dense>self._new_c([zero,one], parent)
        cdef SkewPolynomial_finite_field_dense v
        if (exp % 2 == 1):
            v = <SkewPolynomial_finite_field_dense>self._new_c([zero,one], parent)
            if self.degree() == 1:
                v._inplace_rrem(self)
        else:
            v = <SkewPolynomial_finite_field_dense>self._new_c([one], parent)
        exp = exp >> 1
        n = 1
        while exp != 0:
            powx._inplace_lmul(powx.conjugate(n))
            powx._inplace_rrem(self)
            n = n << 1
            if (exp % 2 == 1):
                v = v.conjugate(n)
                v._inplace_rmul(powx)
                v._inplace_rrem(self)
            exp = exp >> 1
        l = v.list()
        for i from 0 <= i < len(l):
            phir.set_unsafe(i,col,l[i])
        for j from col+1 <= j < d:
            v <<= 1
            v = v.conjugate(1)
            v._inplace_rrem(self)
            for i from 0 <= i <= v.degree():
                phir.set_unsafe(i,j,v.__coeffs[i])
        return phir

    def smurf(self):
        return self._matphir_c()


    cdef Matrix_dense _matmul_c(self):
        r"""
        Return the matrix of the multiplication by self on
        `K[X,\sigma]` considered as a free module over `K[X^r]`
        (here `r` is the order of `\sigma`).

        .. WARNING::

            Does not work if self is not monic.
        """
        cdef Py_ssize_t i, j, deb, k, r = self.parent()._order
        cdef Py_ssize_t d = self.degree ()
        cdef Ring base_ring = <Ring?>self.parent().base_ring()
        cdef RingElement minusone = <RingElement?>base_ring(-1)
        cdef RingElement zero = <RingElement?>base_ring(0)
        cdef Polk = PolynomialRing (base_ring, 'xr')
        cdef Matrix_dense M = <Matrix_dense?>zero_matrix(Polk,r,r)
        cdef list l = self.list()
        for j from 0 <= j < r:
            for i from 0 <= i < r:
                if i < j:
                    pol = [ zero ]
                    deb = i-j+r
                else:
                    pol = [ ]
                    deb = i-j
                for k from deb <= k <= d by r:
                    pol.append(l[k])
                M.set_unsafe(i,j,Polk(pol))
            for i from 0 <= i <= d:
                l[i] = self._parent.twist_map()(l[i])
        return M


    def reduced_norm(self):
        r"""
        Return the reduced norm of this skew polynomial.

        .. NOTE::

            The result is cached.

        ALGORITHM:

        If `r` (= the order of the twist map) is small compared
        to `d` (= the degree of this skew polynomial), the reduced
        norm is computed as the determinant of the multiplication
        by `P` (= this skew polynomial) acting on `K[X,\sigma]`
        (= the underlying skew ring) viewed as a free module of
        rank `r` over `K[X^r]`.

        Otherwise, the reduced norm is computed as the characteristic
        polynomial (considered as a polynomial of the variable `X^r`)
        of the left multiplication by `X` on the quotient
        `K[X,\sigma] / K[X,\sigma]*P` (which is a `K`-vector space
        of dimension `d`).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = S.random_element(degree=3,monic=True); a
            x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: N = a.reduced_norm(); N
            (x^3)^3 + 4*(x^3)^2 + 4

        Note that the parent of `N` is the center of the `S`
        (and not `S` itself)::

            sage: N.parent()
            Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5
            sage: N.parent() == S.center()
            True

        In any case, coercion works fine::

            sage: S(N)
            x^9 + 4*x^6 + 4
            sage: N + a
            x^9 + 4*x^6 + x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 1

        We check that `N` is a multiple of `a`::

            sage: S(N).is_divisible_by(a)
            True
            sage: S(N).is_divisible_by(a,side=Left)
            True

        .. NOTE::

            We really need to coerce first `N` into `S`. Otherwise an
            error occurs::

                 sage: N.is_divisible_by(a)
                 Traceback (most recent call last):
                 ...
                 AttributeError: 'sage.rings.polynomial.skew_polynomial_element.CenterSkewPolynomial_generic_dense' object has no attribute 'is_divisible_by'

        We check that the reduced norm is a multiplicative map::

            sage: a = S.random_element(degree=5)
            sage: b = S.random_element(degree=7)
            sage: a.reduced_norm() * b.reduced_norm() == (a*b).reduced_norm()
            True
        """
#        norm = self.norm()
        if self._norm is None:
#        if norm is None:
            center = self.parent().center()
            if self.is_zero():
#                norm = center(0)
                self._norm = center(0)
            else:
                section = center._embed_basering.section()
                exp = (self.parent().base_ring().cardinality() - 1) / (center.base_ring().cardinality() - 1)
                order = self.parent()._order
                lc = section(self.leading_coefficient()**exp)
                if order < self.degree():
                    M = self._matmul_c()
#                    norm = center([ lc*section(x) for x in M.determinant().monic().list() ])
                    self._norm = center([ lc*section(x) for x in M.determinant().monic().list() ])
                else:
                    charpoly = self._matphir_c().characteristic_polynomial()
#                    norm = center([ lc*section(x) for x in charpoly.list() ])
                    self._norm = center([ lc*section(x) for x in charpoly.list() ])
        return self._norm
#        return norm


    def reduced_norm_factor(self):
        """
        Return the reduced norm of this polynomial
        factorized in the centre.

        EXAMPLES:

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: a = (x^2 + 1) * (x+3)
            sage: a.reduced_norm_factor()
            ((x^3) + 3) * ((x^3) + 2)^2
        """
        if self._norm_factor is None:
            self._norm_factor = self.reduced_norm().factor()
        return self._norm_factor


    def is_central(self):
        """
        Return True if this skew polynomial lies in the center.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: x.is_central()
            False
            sage: (t*x^3).is_central()
            False
            sage: (x^6 + x^3).is_central()
            True
        """
        center = self.parent().center()
        try:
            center(self)
            return True
        except ValueError:
            return False


    def bound(self):
        """
        Return a bound of this skew polynomial (i.e. a multiple
        of this skew polynomial lying in the center).

        .. NOTE::

            Since `b` is central, it divides a skew polynomial
            on the left iff it divides it on the right

        ALGORITHM:

        #. Sage first checks whether ``self`` is itself in the
           center. It if is, it returns ``self``

        #. If an optimal bound was previously computed and
           cached, Sage returns it

        #. Otherwise, Sage returns the reduced norm of ``self``

        As a consequence, the output of this function may depend
        on previous computations (an example is given below).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center()

            sage: a = x^2 + (4*t + 2)*x + 4*t^2 + 3
            sage: b = a.bound(); b
            (x^3)^2 + (x^3) + 4

        Note that the parent of `b` is the center of the `S`
        (and not `S` itself)::

            sage: b.parent()
            Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5
            sage: b.parent() == Z
            True

        We check that `b` is divisible by `a`::

            sage: S(b).is_divisible_by(a)
            True
            sage: S(b).is_divisible_by(a,side=Left)
            True

        Actually, `b` is the reduced norm of `a`::

            sage: b == a.reduced_norm()
            True

        Now, we compute the optimal bound of `a` and see that
        it affects the behaviour of ``bound()``::

            sage: a.optimal_bound()
            (x^3) + 3
            sage: a.bound()
            (x^3) + 3

        We finally check that if `a` is a central skew polynomial,
        then ``a.bound()`` returns simply `a`::

            sage: a = S(Z.random_element(degree=4)); a
            2*x^12 + x^9 + 2*x^3
            sage: b = a.bound(); b
            2*(x^3)^4 + (x^3)^3 + 2*(x^3)
            sage: a == b
            True
        """
        center = self.parent().center()
        try:
            return center(self)
        except ValueError:
            pass
        if not self._optbound is None:
            return center(self._optbound)
        return self.reduced_norm()


    def optimal_bound(self):
        """
        Return the optimal bound of this skew polynomial (i.e.
        the monic multiple of this skew polynomial of minimal
        degree lying in the center).

        .. NOTE::

            The result is cached.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center()

            sage: a = x^2 + (4*t + 2)*x + 4*t^2 + 3
            sage: b = a.optimal_bound(); b
            (x^3) + 3

        Note that the parent of `b` is the center of the `S`
        (and not `S` itself)::

            sage: b.parent()
            Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5
            sage: b.parent() == Z
            True

        We check that `b` is divisible by `a`::

            sage: S(b).is_divisible_by(a)
            True
            sage: S(b).is_divisible_by(a,side=Left)
            True
        """
        center = self.parent().center()
        if self._optbound is None:
            try:
                self._optbound = center(self).monic()
            except ValueError:
                bound = self._matphir_c().minimal_polynomial()
                section = center._embed_basering.section()
                self._optbound = [ section(x) for x in bound.list() ]
        return center(self._optbound)


    def is_irreducible(self):
        """
        Return true if this skew polynomial is irreducible.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: a = x^2 + t*x + 1
            sage: a.is_irreducible()
            False
            sage: a.factor()
            (x + 4*t^2 + 4*t + 1) * (x + 3*t + 2)

            sage: a = x^2 + t*x + t + 1
            sage: a.is_irreducible()
            True
            sage: a.factor()
            x^2 + t*x + t + 1

        Skew polynomials of degree `1` are of course irreducible::

            sage: a = x + t
            sage: a.is_irreducible()
            True

        A random irreducible skew polynomial is irreducible::

            sage: a = S.random_irreducible(degree=4,monic=True); a   # random
            x^4 + (t + 1)*x^3 + (3*t^2 + 2*t + 3)*x^2 + 3*t*x + 3*t
            sage: a.is_irreducible()
            True

        By convention, constant skew polynomials are not irreducible::

            sage: S(1).is_irreducible()
            False
            sage: S(0).is_irreducible()
            False
        """
        return self.reduced_norm().is_irreducible()


    def type(self,N):
        """
        INPUT:

        -  ``N`` -- an irreducible polynomial in the
           center of the underlying skew polynomial ring

        OUTPUT:

        The `N`-type of this skew polynomial

        .. NOTE::

            The result is cached.

        DEFINITION:

        The `N`-type of a skew polynomial `a` is the Partition
        `(t_0, t_1, t_2, ...)` defined by

        .. MATH::

            t_0 + \cdots + t_i = \frac{\deg gcd(a,N^i)}{\deg N}

        where `\deg N` is the degree of `N` considered as an
        element in the center.

        This notion has an important mathematic interest because
        it corresponds to the Jordan type of the `N`-typical part
        of the associated Galois representation.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center(); x3 = Z.gen()

            sage: a = x^4 + x^3 + (4*t^2 + 4)*x^2 + (t^2 + 2)*x + 4*t^2
            sage: N = x3^2 + x3 + 1
            sage: a.type(N)
            [1]
            sage: N = x3 + 1
            sage: a.type(N)
            [2]

            sage: a = x^3 + (3*t^2 + 1)*x^2 + (3*t^2 + t + 1)*x + t + 1
            sage: N = x3 + 1
            sage: a.type(N)
            [2, 1]

        If `N` does not divide the reduced map of `a`, the type
        is empty::

            sage: N = x3 + 2
            sage: a.type(N)
            []

        If `a = N`, the type is just `[r]` where `r` is the order
        of the twist map ``Frob``::

            sage: N = x3^2 + x3 + 1
            sage: S(N).type(N)
            [3]

        `N` must be irreducible::

            sage: N = (x3 + 1) * (x3 + 2)
            sage: a.type(N)
            Traceback (most recent call last):
            ...
            ValueError: N is not irreducible
        """
        try:
            return self._types[N]
        except (KeyError, TypeError):
            if not N.is_irreducible():
                raise ValueError("N is not irreducible")
            skew_ring = self._parent
            if self._norm_factor is None:
                m = -1
            else:
                i = [ n for n,_ in self._norm_factor ].index(N)
                m = self._norm_factor[i][1]
            NS = skew_ring(N)
            type = [ ]
            degN = N.degree()
            while True:
                d = self.gcd(NS)
                deg = d.degree()/degN
                if deg == 0:
                    break
                if m >= 0:
                    if deg == 1:
                        type += m * [1]
                        break
                    m -= deg
                self = self // d
                type.append(deg)
            type = Partition(type)
            if self._types is None:
                self._types = { N: type }
            else:
                self._types[N] = type
            return type


    # Finding divisors
    # ----------------

    cdef SkewPolynomial_finite_field_dense _rdivisor_c(P, CenterSkewPolynomial_generic_dense N):
        """
        cython procedure computing an irreducible monic right divisor
        of `P` whose reduced norm is `N`

        .. WARNING::

            `N` needs to be an irreducible factor of the
            reduced norm of `P`. This function does not check
            this (and his behaviour is not defined if the
            require property doesn't hold).
        """
        cdef skew_ring = P._parent
        cdef Py_ssize_t d = N.degree()
        cdef Py_ssize_t e = P.degree()/d
        cdef SkewPolynomial_finite_field_dense D
        if e == 1:
            D = <SkewPolynomial_finite_field_dense>P._new_c(list(P.__coeffs),skew_ring)
            D._inplace_rmonic()
            return D

        E = N.parent().base_ring().extension(N,name='xr')
        PE = PolynomialRing(E,name='T')
        cdef Integer exp
        if skew_ring.characteristic() != 2:
            exp = Integer((E.cardinality()-1)/2)
        cdef SkewPolynomial_finite_field_dense NS = <SkewPolynomial_finite_field_dense>skew_ring(N)
        cdef SkewPolynomial_finite_field_dense Q = <SkewPolynomial_finite_field_dense>(NS // P)
        cdef SkewPolynomial_finite_field_dense R, X
        cdef Matrix_dense M = <Matrix_dense?>MatrixSpace(E,e,e)(0)
        cdef Matrix_dense V = <Matrix_dense?>MatrixSpace(E,e,1)(0)
        cdef Matrix_dense W
        cdef Py_ssize_t i, j, t, r = skew_ring._order
        cdef Polynomial dd, xx, yy, zz

        while 1:
            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((e*r-1,e*r-1))
            R._inplace_lmul(Q)
            X = <SkewPolynomial_finite_field_dense>Q._new_c(Q.__coeffs[:],Q._parent)
            for j from 0 <= j < e:
                for i from 0 <= i < e:
                    M.set_unsafe(i, j, E([skew_ring._retraction(X[t*r+i]) for t in range(d)]))
                X._inplace_lmul(R)
                X._inplace_rrem(NS)
            for i from 0 <= i < e:
                V.set_unsafe(i, 0, E([skew_ring._retraction(X[t*r+i]) for t in range(d)]))
            W = M._solve_right_nonsingular_square(V)
            if M*W != V:
                skew_ring._new_retraction_map()
                continue
            xx = PE(W.list()+[E(-1)])
            if skew_ring.characteristic() == 2:
                yy = PE.gen()
                zz = PE.gen()
                for i from 1 <= i < d:
                    zz = (zz*zz) % xx
                    yy += zz
                dd = xx.gcd(yy)
                if dd.degree() != 1: continue
            else:
                yy = PE.gen().__pow__(exp,xx) - 1
                dd = xx.gcd(yy)
                if dd.degree() != 1:
                    yy += 2
                    dd = xx.gcd(yy)
                    if dd.degree() != 1: continue
            D = P._rgcd(R + skew_ring.center()((dd[0]/dd[1]).list()))
            if D.degree() == 0:
                continue
            D._inplace_rmonic()
            D._init_cache()
            return D


    def irreducible_divisor(self,side=Right,distribution=None):
        """
        INPUT:

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        -  ``distribution`` -- None (default) or ``uniform``

           - None: no particular specification

           - ``uniform``: the returned irreducible divisor is
             uniformly distributed

        .. NOTE::

            ``uniform`` is a little bit slower.

        OUTPUT:

        -  an irreducible monic ``side`` divisor of ``self``

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^6 + 3*t*x^5 + (3*t + 1)*x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + 2)*x + 4*t^2 + 3*t + 3

            sage: dr = a.irreducible_divisor(); dr  # random
            x^3 + (2*t^2 + t + 4)*x^2 + (4*t + 1)*x + 4*t^2 + t + 1
            sage: a.is_divisible_by(dr)
            True

            sage: dl = a.irreducible_divisor(side=Left); dl  # random
            x^3 + (2*t^2 + t + 1)*x^2 + (4*t^2 + 3*t + 3)*x + 4*t^2 + 2*t + 1
            sage: a.is_divisible_by(dl,side=Left)
            True

        Right divisors are cached. Hence, if we ask again for a
        right divisor, we will get the same answer::

            sage: a.irreducible_divisor()  # random
            x^3 + (2*t^2 + t + 4)*x^2 + (4*t + 1)*x + 4*t^2 + t + 1

        However the algorithm is probabilistic. Hence, if we first
        reinitialiaze `a`, we may get a different answer::

            sage: a = x^6 + 3*t*x^5 + (3*t + 1)*x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + 2)*x + 4*t^2 + 3*t + 3
            sage: a.irreducible_divisor()  # random
            x^3 + (t^2 + 3*t + 4)*x^2 + (t + 2)*x + 4*t^2 + t + 1

        We can also generate uniformly distributed irreducible monic
        divisors as follows::

            sage: a.irreducible_divisor(distribution="uniform")  # random
            x^3 + (4*t + 2)*x^2 + (2*t^2 + 2*t + 2)*x + 2*t^2 + 2
            sage: a.irreducible_divisor(distribution="uniform")  # random
            x^3 + (t^2 + 2)*x^2 + (3*t^2 + 1)*x + 4*t^2 + 2*t
            sage: a.irreducible_divisor(distribution="uniform")  # random
            x^3 + x^2 + (4*t^2 + 2*t + 4)*x + t^2 + 3

        By convention, the zero skew polynomial has no irreducible
        divisor:

            sage: S(0).irreducible_divisor()
            Traceback (most recent call last):
            ...
            ValueError: 0 has no irreducible divisor
        """
        if self.is_zero():
            raise ValueError("0 has no irreducible divisor")
        if not (distribution is None or distribution == "uniform"):
            raise ValueError("distribution must be None or 'uniform'")
        if distribution == "uniform":
            skew_ring = self._parent
            center = skew_ring.center()
            cardcenter = center.base_ring().cardinality()
            gencenter = center.gen()
            count = [ ]
            total = 0
            F = self.reduced_norm_factor()
            for n,_ in F:
                if n == gencenter:
                    total += 1
                else:
                    degn = n.degree()
                    P = self.gcd(skew_ring(n))
                    m = P.degree()/degn
                    cardL = cardcenter**degn
                    total += (cardL**m - 1)/(cardL - 1)
                count.append(total)
            if total == 0:
                raise ValueError("No irreducible divisor having given reduced norm")
            random = ZZ.random_element(total)
            for i in range(len(F)):
                if random < count[i]:
                    N = F[i][0]
                    break
        else:
            N = self.reduced_norm_factor()[0][0]
        return self.irreducible_divisor_with_norm(N,side=side,distribution=distribution)


    def irreducible_divisor_with_norm(self,N,side=Right,distribution=None): # Ajouter side
        """
        INPUT:

        -  ``N`` -- an irreducible polynomial in the center
           of the underlying skew polynomial ring

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        -  ``distribution`` -- None (default) or ``uniform``

           - None: no particular specification

           - ``uniform``: the returned irreducible divisor is
             uniformly distributed

        .. NOTE::

            ``uniform`` is a little bit slower.

        OUTPUT:

        -  an irreducible monic ``side`` divisor of ``self``
           whose reduced norm is similar to `N` (i.e. `N` times
           a unit).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center(); x3 = Z.gen()
            sage: a = x^6 + 3*x^3 + 2

            sage: d1 = a.irreducible_divisor_with_norm(x3+1); d1   # random
            x + t^2 + 3*t
            sage: a.is_divisible_by(d1)
            True
            sage: d1.reduced_norm()
            (x^3) + 1

            sage: d2 = a.irreducible_divisor_with_norm(x3+2); d2   # random
            x + 2*t^2 + 3*t + 2
            sage: a.is_divisible_by(d2)
            True
            sage: d2.reduced_norm()
            (x^3) + 2

            sage: d3 = a.irreducible_divisor_with_norm(x3+3)
            Traceback (most recent call last):
            ...
            ValueError: No irreducible divisor having given reduced norm

        We can also generate uniformly distributed irreducible monic
        divisors as follows::

            sage: a.irreducible_divisor_with_norm(x3+1,distribution="uniform")   # random
            x + 3*t^2 + 3*t + 1
            sage: a.irreducible_divisor_with_norm(x3+1,distribution="uniform")   # random
            x + 1
            sage: a.irreducible_divisor_with_norm(x3+1,distribution="uniform")   # random
            x + 2*t^2 + 4*t
        """
        cdef SkewPolynomial_finite_field_dense cP1
        cdef CenterSkewPolynomial_generic_dense cN
        if self.is_zero():
            raise "No irreducible divisor having given reduced norm"
        skew_ring = self._parent
        center = skew_ring.center()
        try:
            N = center(N)
        except TypeError:
            raise TypeError("N must be a polynomial in the center")
        cardcenter = center.base_ring().cardinality()
        gencenter = center.gen()

        if N == gencenter:
            if self[0] == 0:
                return skew_ring.gen()
            else:
                raise ValueError("No irreducible divisor having given reduced norm")

        D = None
        try:
            D = self._rdivisors[N]
        except (KeyError, TypeError):
            if N.is_irreducible():
                cP1 = <SkewPolynomial_finite_field_dense>self._rgcd(self._parent(N))
                cN = <CenterSkewPolynomial_generic_dense>N
                if cP1.degree() > 0:
                    D = cP1._rdivisor_c(cN)
            if self._rdivisors is None:
                self._rdivisors = { N: D }
            else:
                self._rdivisors[N] = D
            distribution = ""
        if D is None:
            raise ValueError("No irreducible divisor having given reduced norm")

        NS = self._parent(N)
        degN = N.degree()
        if side is Right:
            if distribution == "uniform":
                P1 = self._rgcd(NS)
                if P1.degree() != degN:
                    Q1 = NS // P1
                    deg = P1.degree()-1
                    while True:
                        R = Q1*skew_ring.random_element((deg,deg))
                        if P1.gcd(R) == 1:
                            break
                    D = P1.gcd(D*R)
            return D
        else:
            deg = NS.degree()-1
            P1 = self.lgcd(NS)
            while True:
                if distribution == "uniform":
                    while True:
                        R = skew_ring.random_element((deg,deg))
                        if NS.gcd(R) == 1:
                            break
                    D = NS.gcd(D*R)
                Dp = NS // D
                LDp = P1.gcd(Dp)
                LD = P1 // LDp
                if LD.degree() == degN:
                    return LD
                distribution = "uniform"


    def irreducible_divisors(self,side=Right):
        """
        INPUT:

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        An iterator over all irreducible monic ``side`` divisors
        of this skew polynomial

        EXAMPLES:

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + 2*t*x^3 + 3*t^2*x^2 + (t^2 + t + 1)*x + 4*t + 3
            sage: iter = a.irreducible_divisors(); iter
            <generator object at 0x...>
            sage: iter.next()   # random
            x + 2*t^2 + 4*t + 4
            sage: iter.next()   # random
            x + 3*t^2 + 4*t + 1

        We can use this function to build the list of all monic
        irreducible divisors of `a`::

            sage: rightdiv = [ d for d in a.irreducible_divisors() ]
            sage: leftdiv = [ d for d in a.irreducible_divisors(side=Left) ]

        We do some checks::

            sage: len(rightdiv) == a.count_irreducible_divisors()
            True
            sage: len(rightdiv) == len(Set(rightdiv))  # check no duplicates
            True
            sage: for d in rightdiv:
            ...       if not a.is_divisible_by(d):
            ...           print "Found %s which is not a right divisor" % d
            ...       elif not d.is_irreducible():
            ...           print "Found %s which is not irreducible" % d

            sage: len(leftdiv) == a.count_irreducible_divisors(side=Left)
            True
            sage: len(leftdiv) == len(Set(leftdiv))  # check no duplicates
            True
            sage: for d in leftdiv:
            ...       if not a.is_divisible_by(d,side=Left):
            ...           print "Found %s which is not a left divisor" % d
            ...       elif not d.is_irreducible():
            ...           print "Found %s which is not irreducible" % d

        Note that left divisors and right divisors differ::

            sage: Set(rightdiv) == Set(leftdiv)
            False

        Note that the algorithm is probabilistic. As a consequence, if we
        build again the list of right monic irreducible divisors of `a`, we
        may get a different ordering::

            sage: rightdiv2 = [ d for d in a.irreducible_divisors() ]
            sage: rightdiv == rightdiv2
            False
            sage: Set(rightdiv) == Set(rightdiv2)
            True
        """
        return self._irreducible_divisors(side)


    def _irreducible_divisors(self,side): # prendre side en compte
        """
        Return an iterator over all irreducible monic
        divisors of this skew polynomial.

        Do not use this function. Use instead
        ``self.irreducible_divisors()``.
        """
        if self.is_zero():
            return
        skew_ring = self._parent
        center = skew_ring.center()
        kfixed = center.base_ring()
        F = self.reduced_norm_factor()
        oppside = side.opposite()
        for N,_ in F:
            if N == center.gen():
                yield skew_ring.gen()
                continue
            degN = N.degree()
            NS = skew_ring(N)
            P = self.gcd(NS,side=side)
            m = P.degree()/degN
            if m == 1:
                yield P
                continue
            degrandom = P.degree() - 1
            Q,_ = NS.quo_rem(P,side=side)
            P1 = self.irreducible_divisor_with_norm(N,side=side)
            Q1,_ = P.quo_rem(P1,side=side)
            while True:
                R = skew_ring.random_element((degrandom,degrandom))
                if side is Right:
                    g = (R*Q).rem(P,side=Left)
                else:
                    g = (Q*R).rem(P)
                if g.gcd(P,side=oppside) != 1: continue
                L = Q1
                V = L
                for i in range(1,m):
                    if side is Right:
                        L = (g*L).gcd(P,side=Left)
                    else:
                        L = (L*g).gcd(P)
                    V = V.gcd(L,side=oppside)
                if V == 1: break
            rng = xmrange_iter([kfixed]*degN,center)
            for i in range(m):
                for pol in xmrange_iter([rng]*i):
                    f = skew_ring(1)
                    for j in range(i):
                        coeff = pol.pop()
                        f = (g*f+coeff).rem(P,side=oppside)
                    if side is Right:
                        d = (f*Q1).gcd(P,side=Left)
                    else:
                        d = (Q1*f).gcd(P)
                    d,_ = P.quo_rem(d,side=oppside)
                    yield d


    def count_irreducible_divisors(self,side=Right):
        """
        INPUT:

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        The number of irreducible monic ``side`` divisors of
        this skew polynomial.

        .. NOTE::

            Actually, one can prove that there are always as
            many left irreducible monic divisors as right
            irreducible monic divisors.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

        We illustrate that a skew polynomial may have a number of irreducible
        divisors greater than its degree.

            sage: a = x^4 + (4*t + 3)*x^3 + t^2*x^2 + (4*t^2 + 3*t)*x + 3*t
            sage: a.count_irreducible_divisors()
            12
            sage: a.count_irreducible_divisors(side=Left)
            12

        We illustrate that an irreducible polynomial in the center have
        in general a lot of irreducible divisors in the skew polynomial
        ring::

            sage: Z = S.center(); x3 = Z.gen()
            sage: N = x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3; N
            (x^3)^5 + 4*(x^3)^4 + 4*(x^3)^2 + 4*(x^3) + 3
            sage: N.is_irreducible()
            True
            sage: S(N).count_irreducible_divisors()
            9768751
        """
        if self.is_zero():
            return 0
        skew_ring = self.parent()
        cardcenter = skew_ring.center().base_ring().cardinality()
        gencenter = skew_ring.center().gen()
        F = self.reduced_norm_factor()
        val = self.valuation()
        self >>= val
        count = 0
        if val > 0:
            count = 1
        for N,_ in F:
            if N == gencenter:
                continue
            degN = N.degree()
            P = self.gcd(skew_ring(N), side=side)
            m = P.degree()/degN
            cardL = cardcenter**degN
            count += (cardL**m - 1)/(cardL - 1)
        return count


    # Finding factorizations
    # ----------------------

    cdef _factor_c(self):
        """
        Compute a factorization of ``self``
        """
        cdef skew_ring = self._parent
        cdef Py_ssize_t degQ, degrandom, m, mP, i
        cdef CenterSkewPolynomial_generic_dense N
        cdef SkewPolynomial_finite_field_dense poly = <SkewPolynomial_finite_field_dense>self.rmonic()
        cdef val = poly._val_inplace_unit()
        if val == -1:
            return Factorization([], sort=False, unit=skew_ring.zero())
#            return Factorization([], sort=False, unit=skew_ring.zero_element())
        cdef list factors = [ (skew_ring.gen(), val) ]
        cdef SkewPolynomial_finite_field_dense P, Q, P1, NS, g, right, Pn
        cdef SkewPolynomial_finite_field_dense right2 = skew_ring(1) << val
        cdef RingElement unit = <RingElement>self.leading_coefficient()
        cdef Polynomial gencenter = skew_ring.center().gen()
        cdef Py_ssize_t p = skew_ring.characteristic()
        cdef F = self.reduced_norm_factor()

        for N,m in F:
            if N == gencenter:
                continue
            degN = N.degree()
            if poly.degree() == degN:
                factors.append((poly,1))
                break
            NS = <SkewPolynomial_finite_field_dense>skew_ring(N)
            P1 = None
            while 1:
                P = <SkewPolynomial_finite_field_dense>poly._rgcd(NS)
                P._inplace_rmonic()
                mP = P.degree() / degN
                if mP == 0: break
                if mP == 1:
                    factors.append((P,1))
                    poly._inplace_rfloordiv(P)
                    for i from 1 <= i < m:
                        if poly.degree() == degN:
                            factors.append((poly,1))
                            break
                        P = poly._rgcd(NS)
                        P._inplace_rmonic()
                        factors.append((P,1))
                        poly._inplace_rfloordiv(P)
                    break
                if P1 is None:
                    P1 = P._rdivisor_c(N)
                Q = <SkewPolynomial_finite_field_dense>NS._new_c(NS.__coeffs[:], NS._parent)
                Q._inplace_rfloordiv(P)
                Q._inplace_lmul(P1)
                factors.append((P1,1))
                right = <SkewPolynomial_finite_field_dense>P1._new_c(P1.__coeffs[:], P1._parent)
                m -= (mP-1)
                degrandom = P.degree()
                while mP > 2:
                    while 1:
                        g = <SkewPolynomial_finite_field_dense>skew_ring.random_element((degrandom,degrandom))
                        g._inplace_lmul(Q)
                        g._inplace_rgcd(P)
                        Pn = right._coeff_llcm(g)
                        if len(Pn.__coeffs)-1 == degN: break
                    Pn._inplace_rmonic()
                    factors.append((Pn,1))
                    right._inplace_lmul(Pn)
                    degrandom -= degN
                    mP -= 1
                poly._inplace_rfloordiv(right)
                P1,_ = P.rquo_rem(right)
        factors.reverse()
        return Factorization(factors, sort=False, unit=unit)


    cdef _factor_uniform_c(self):
        """
        Compute a uniformly distrbuted factorization of ``self``
        """
        skew_ring = self._parent
        cdef Integer cardE, cardcenter = skew_ring.center().base_ring().cardinality()
        cdef CenterSkewPolynomial_generic_dense gencenter = <CenterSkewPolynomial_generic_dense>skew_ring.center().gen()
        cdef SkewPolynomial_finite_field_dense gen = <SkewPolynomial_finite_field_dense>skew_ring.gen()

        cdef list factorsN = [ ]
        cdef dict dict_divisor = { }
        cdef dict dict_type = { }
        cdef dict dict_right = { }
        cdef CenterSkewPolynomial_generic_dense N
        cdef Py_ssize_t m
        cdef list type

        for N,m in self.reduced_norm_factor():
            factorsN += m * [N]
            if N == gencenter: continue
            type = list(self.type(N))
            dict_type[N] = type
            if type[0] > 1:
                dict_divisor[N] = self.irreducible_divisor_with_norm(N)
                dict_right[N] = skew_ring(1)
        cdef list indices = list(Permutations(len(factorsN)).random_element())

        cdef RingElement unit = self.leading_coefficient()
        cdef SkewPolynomial_finite_field_dense left = self._new_c(self.__coeffs[:],skew_ring)
        left._inplace_rmonic()
        cdef SkewPolynomial_finite_field_dense right = <SkewPolynomial_finite_field_dense>skew_ring(1)
        cdef SkewPolynomial_finite_field_dense L, R
        cdef SkewPolynomial_finite_field_dense NS, P, Q, D, D1, D2, d
        cdef list factors = [ ]
        cdef list maxtype
        cdef Py_ssize_t i, j, degN, deg
        cdef count, maxcount

        for i in indices:
            N = factorsN[i-1]
            if N == gencenter:
                D1 = gen
            else:
                type = dict_type[N]
                NS = skew_ring(N)
                P = left.gcd(NS)
                if type[0] == 1:
                    D1 = P
                else:
                    R = right._new_c(right.__coeffs[:],skew_ring)
                    R._inplace_rfloordiv(dict_right[N])
                    D = R._coeff_llcm(dict_divisor[N])
                    maxtype = list(type)
                    maxtype[-1] -= 1
                    degN = N.degree()
                    cardE = cardcenter ** degN
                    maxcount = q_jordan(Partition(maxtype),cardE)
                    Q = NS // P
                    deg = P.degree()-1
                    while 1:
                        while 1:
                            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((deg,deg))
                            R._inplace_lmul(Q)
                            if P._rgcd(R).degree() == 0:
                                break
                        D1 = P._rgcd(D*R)
                        D1._inplace_rmonic()

                        L = left._new_c(list(left.__coeffs),skew_ring)
                        L._inplace_rfloordiv(D1)
                        degN = N.degree()
                        for j in range(len(type)):
                            if type[j] == 1:
                                newtype = type[:-1]
                                break
                            d = L._rgcd(NS)
                            d._inplace_rmonic()
                            deg = d.degree() / degN
                            if deg < type[j]:
                                newtype = type[:]
                                newtype[j] = deg
                                break
                            L._inplace_rfloordiv(d)
                        count = q_jordan(Partition(newtype),cardE)
                        if ZZ.random_element(maxcount) < count:
                            break
                    dict_type[N] = newtype

                    D2 = D._new_c(list(D.__coeffs),skew_ring)
                    D2._inplace_rmonic()
                    while D2 == D1:
                        while 1:
                            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((deg,deg))
                            R._inplace_lmul(Q)
                            if P._rgcd(R).degree() == 0:
                                break
                        D2 = P._rgcd(D*R)
                        D2._inplace_rmonic()
                    dict_divisor[N] = D1._coeff_llcm(D2)
            factors.append((D1,1))
            left._inplace_rfloordiv(D1)
            right._inplace_lmul(D1)
            dict_right[N] = right._new_c(list(right.__coeffs),skew_ring)

        factors.reverse()
        return Factorization(factors,sort=False,unit=unit)


    def factor(self,distribution=None):
        """
        Return a factorization of this skew polynomial.

        INPUT:

        -  ``distribution`` -- None (default) or ``uniform``

           - None: no particular specification

           - ``uniform``: the returned factorization is uniformly
             distributed among all possible factorizations

        .. NOTE::

            ``uniform`` is a little bit slower.

        OUTPUT:

        -  ``Factorization`` -- a factorization of self as a
           product of a unit and a product of irreducible monic
           factors

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^3 + (t^2 + 4*t + 2)*x^2 + (3*t + 3)*x + t^2 + 1
            sage: F = a.factor(); F  # random
            (x + t^2 + 4) * (x + t + 3) * (x + t)
            sage: F.value() == a
            True

        The result of the factorization is cached. Hence, if we try
        again to factor `a`, we will get the same answer::

            sage: a.factor()  # random
            (x + t^2 + 4) * (x + t + 3) * (x + t)

        However, the algorithm is probabilistic. Hence if we first
        reinitialiaze `a`, we may get a different answer::

            sage: a = x^3 + (t^2 + 4*t + 2)*x^2 + (3*t + 3)*x + t^2 + 1
            sage: F = a.factor(); F   # random
            (x + t^2 + t + 2) * (x + 2*t^2 + t + 4) * (x + t)
            sage: F.value() == a
            True

        There is no guarantee on the distribution of the factorizations
        we get that way. (For this particular `a` for example, we get the
        uniform distribution on the subset of all factorizations ending
        by the factor `x + t`.)

        If we rather want uniform distribution among all factorizations,
        we need to specify it as follows::

            sage: a.factor(distribution="uniform")   # random
            (x + t^2 + 4) * (x + t) * (x + t + 3)
            sage: a.factor(distribution="uniform")   # random
            (x + 2*t^2) * (x + t^2 + t + 1) * (x + t^2 + t + 2)
            sage: a.factor(distribution="uniform")   # random
            (x + 2*t^2 + 3*t) * (x + 4*t + 2) * (x + 2*t + 2)

        By convention, the zero skew polynomial has no factorization:

            sage: S(0).factor()
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined
        """
        if not (distribution is None or distribution == "uniform"):
            raise ValueError("distribution must be None or 'uniform'")
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        sig_on()
        if distribution is None:
            if self._factorization is None:
                self._factorization = self._factor_c()
            F = self._factorization
        else:
            F = self._factor_uniform_c()
            if self._factorization is None:
                self._factorization = F
        sig_off()
        return F


    def count_factorizations(self):
        """
        Return the number of factorizations (as a product of a
        unit and a product of irreducible monic factors) of this
        skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + (4*t + 3)*x^3 + t^2*x^2 + (4*t^2 + 3*t)*x + 3*t
            sage: a.count_factorizations()
            216

        We illustrate that an irreducible polynomial in the center have
        in general a lot of distinct factorizations in the skew polynomial
        ring::

            sage: Z = S.center(); x3 = Z.gen()
            sage: N = x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3; N
            (x^3)^5 + 4*(x^3)^4 + 4*(x^3)^2 + 4*(x^3) + 3
            sage: N.is_irreducible()
            True
            sage: S(N).count_factorizations()
            30537115626
        """
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        cardcenter = self._parent.center().base_ring().cardinality()
        gencenter = self._parent.center().gen()
        F = self.reduced_norm_factor()
        summ = 0
        count = 1
        for N,m in F:
            summ += m
            if m == 1: continue
            if N != gencenter:
                count *= q_jordan(self.type(N),cardcenter**N.degree())
            count /= factorial(m)
        return count * factorial(summ)

    def count_factorisations(self):
        """
        Return the number of factorisations (as a product of a
        unit and a product of irreducible monic factors) of this
        skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + (4*t + 3)*x^3 + t^2*x^2 + (4*t^2 + 3*t)*x + 3*t
            sage: a.count_factorisations()
            216

        We illustrate that an irreducible polynomial in the center have
        in general a lot of distinct factorisations in the skew polynomial
        ring::

            sage: Z = S.center(); x3 = Z.gen()
            sage: N = x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3; N
            (x^3)^5 + 4*(x^3)^4 + 4*(x^3)^2 + 4*(x^3) + 3
            sage: N.is_irreducible()
            True
            sage: S(N).count_factorisations()
            30537115626
        """
        return self.count_factorizations()


    # Not optimized (many calls to reduced_norm, reduced_norm_factor,_rdivisor_c, which are slow)
    def factorizations(self):
        """
        Return an iterator over all factorizations (as a product
        of a unit and a product of irreducible monic factors) of
        this skew polynomial.

        .. NOTE::

            The algorithm is probabilistic. As a consequence, if
            we execute two times with the same input we can get
            the list of all factorizations in two differents orders.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^3 + (t^2 + 1)*x^2 + (2*t + 3)*x + t^2 + t + 2
            sage: iter = a.factorizations(); iter
            <generator object at 0x...>
            sage: iter.next()   # random
            (x + 3*t^2 + 4*t) * (x + 2*t^2) * (x + 4*t^2 + 4*t + 2)
            sage: iter.next()   # random
            (x + 3*t^2 + 4*t) * (x + 3*t^2 + 2*t + 2) * (x + 4*t^2 + t + 2)

        We can use this function to build the list of factorizations
        of `a`::

            sage: factorizations = [ F for F in a.factorizations() ]

        We do some checks::

            sage: len(factorizations) == a.count_factorizations()
            True
            sage: len(factorizations) == len(Set(factorizations))  # check no duplicates
            True
            sage: for F in factorizations:
            ...       if F.value() != a:
            ...           print "Found %s which is not a correct factorization" % d
            ...           continue
            ...       for d,_ in F:
            ...           if not d.is_irreducible():
            ...               print "Found %s which is not a correct factorization" % d
        """
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        unit = self.leading_coefficient()
        poly = self.rmonic()
        for factors in self._factorizations_rec():
            yield Factorization(factors,sort=False,unit=unit)
    def _factorizations_rec(self):
        if self.is_irreducible():
            yield [ (self,1) ]
        else:
            for div in self._irreducible_divisors(Right):
                poly = self // div
                # Here, we should update poly._norm, poly._norm_factor, poly._rdivisors
                for factors in poly._factorizations_rec():
                    factors.append((div,1))
                    yield factors


    def factorisations(self):
        """
        Return an iterator over all factorisations (as a product
        of a unit and a product of irreducible monic factors) of
        this skew polynomial.

        .. NOTE::

            The algorithm is probabilistic. As a consequence, if
            we execute two times with the same input we can get
            the list of all factorizations in two differents orders.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^3 + (t^2 + 1)*x^2 + (2*t + 3)*x + t^2 + t + 2
            sage: iter = a.factorisations(); iter
            <generator object at 0x...>
            sage: iter.next()   # random
            (x + 3*t^2 + 4*t) * (x + 2*t^2) * (x + 4*t^2 + 4*t + 2)
            sage: iter.next()   # random
            (x + 3*t^2 + 4*t) * (x + 3*t^2 + 2*t + 2) * (x + 4*t^2 + t + 2)

        We can use this function to build the list of factorizations
        of `a`::

            sage: factorisations = [ F for F in a.factorisations() ]

        We do some checks::

            sage: len(factorisations) == a.count_factorisations()
            True
            sage: len(factorisations) == len(Set(factorisations))  # check no duplicates
            True
            sage: for F in factorisations:
            ...       if F.value() != a:
            ...           print "Found %s which is not a correct factorization" % d
            ...           continue
            ...       for d,_ in F:
            ...           if not d.is_irreducible():
            ...               print "Found %s which is not a correct factorization" % d
        """
        return self.factorizations()


# Karatsuba class
#################

cdef class SkewPolynomial_finite_field_karatsuba:
    """
    A special class implementing Karatsuba multiplication and
    euclidean division on skew polynomial rings over finite fields.

    Only for internal use!

    TESTS::

        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]

    An instance of this class is automatically created when the skew
    ring is created. It is accessible via the key work ``_karatsuba_class``::

        sage: KarClass = S._karatsuba_class; KarClass
        Special class for Karatsuba multiplication and euclidean division on
        Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5

    Each instance of this class is equipped with a ``cutoff`` variable. The
    default value is the maximum between 150 and the order of the twist map
    acting on the finite field::

        sage: KarClass.get_cutoff()
        150

    We can set a new value to this variable as follows::

        sage: KarClass.set_cutoff(27)
        sage: KarClass.get_cutoff()
        27
    """
    def __init__(self,parent,cutoff=0):
        """
        Initialize a new instance of this class.
        """
        self._parent = parent
        self._order = parent._order
        self._zero = parent.base_ring()(0)
        self.set_cutoff(cutoff)
        self._algo_matrix = 0
        self._t = None
        self._T = None
        self._Tinv = None


    def __repr__(self):
        """
        Return a representation of ``self``
        """
        return "Special class for Karatsuba multiplication and euclidean division on\n%s" % self._parent


    def set_cutoff(self,cutoff=0):
        """
        Set a new cutoff for all Karatsuba multiplications
        and euclidean divisions performed by this class.

        INPUT:

        -  ``cutoff`` -- a nonnegative integer or +Infinity
           (default: 0)

        .. NOTE::

            If ``cutoff`` is `0`, the new cutoff is set to the
            default value which is the maximum between 150 and
            the order of the twist map acting on the underlying
            finite field.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: KarClass = S._karatsuba_class

            sage: KarClass.set_cutoff(27)
            sage: KarClass.get_cutoff()
            27

            sage: KarClass.set_cutoff(Infinity)
            sage: KarClass.get_cutoff()
            +Infinity

            sage: KarClass.set_cutoff(0)  # set the default value
            sage: KarClass.get_cutoff()
            150

        The cutoff can't be less than the order of the twist map::

            sage: KarClass.set_cutoff(2)
            Traceback (most recent call last):
            ...
            ValueError: cutoff must be 0 or >= 3
        """
        if cutoff == 0:
            self._cutoff = max(150,self._order)
        else:
            if cutoff < self._order:
                raise ValueError("cutoff must be 0 or >= %s" % self._order)
            if cutoff == Infinity:
                self._cutoff = -1
            else:
                self._cutoff = cutoff


    def get_cutoff(self):
        """
        Return the current cutoff

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: KarClass = S._karatsuba_class
            sage: KarClass.get_cutoff()
            150
            sage: KarClass.set_cutoff(27)
            sage: KarClass.get_cutoff()
            27
        """
        if self._cutoff < 0:
            return Infinity
        else:
            return self._cutoff


    def mul(self,left,right):
        """
        INPUT:

        -  ``left`` -- a skew polynomial in the skew polynomial
           right attached to the class

        -  ``right`` -- a skew polynomial in the skew polynomial
           right attached to the class

        OUTPUT:

        The product left * right (computed by Karatsuba's algorithm)

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: KarClass = S._karatsuba_class

            sage: a = S.random_element(degree=500)
            sage: b = S.random_element(degree=500)
            sage: c = KarClass.mul(a,b)
            sage: c == a*b
            True

        .. NOTE::

            Behind the scene, the operator ``*`` calls actually
            the function ``KarClass.mul()``.
        """
        cdef Py_ssize_t dx = left.degree(), dy = right.degree()
        cdef list x = left.list()
        cdef list y = right.list()
        cdef list res
        if self._cutoff < 0:
            res = self.mul_step(x,y)
        else:
            if dx < dy:
                res = self.mul_iter(y,x,1)
            else:
                res = self.mul_iter(x,y,0)
        return self._parent(res)


    def mul_matrix(self,left,right):
        """
        INPUT:

        -  ``left`` -- an element in the skew polynomial ring
           attached to the class

        -  ``right`` -- an element in the skew polynomial ring
           attached to the class

        OUTPUT:

        The product left * right (computed by Karatsuba-matrix
        algorithm)

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: KarClass = S._karatsuba_class

            sage: a = S.random_element(degree=500)
            sage: b = S.random_element(degree=500)
            sage: c = KarClass.mul_matrix(a,b)
            sage: c == a*b
            True
        """
        cdef Py_ssize_t dx = left.degree(), dy = right.degree()
        cdef list x = left.list()
        cdef list y = right.list()
        cdef list res
        cdef Py_ssize_t save_cutoff = self._cutoff
        cdef Py_ssize_t i, j
        self._cutoff = int(self._order * self._order / 2)
        self._algo_matrix = 1
        if self._t is None:
            self._t = self._parent.base_ring().gen()
        if self._T is None:
            self._T = <Matrix_dense>zero_matrix(self._parent.base_ring(),self._order)
            for i from 0 <= i < self._order:
                map = self._parent.twist_map(-i)
                for j from 0 <= j < self._order:
                    self._T.set_unsafe(i,j,map(self._t**j))
        if self._Tinv is None:
            self._Tinv = <Matrix_dense>self._T.inverse()
        if dx < dy:
            res = self.mul_iter(y,x,1)
        else:
            res = self.mul_iter(x,y,0)
        self._cutoff = save_cutoff
        self._algo_matrix = 0
        return self._parent(res)


    def mul_list(self,x,y):
        """
        INPUT:

        -  ``x`` -- a list of coefficients

        -  ``y`` -- a list of coefficients

        OUTPUT:

        The list of coefficients of the product `a * b`
        where `a` (resp. `b`) is the skew polynomial whose
        coefficients are given by `x` (resp. `y`).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<u> = k['u',Frob]
            sage: KarClass = S._karatsuba_class

            sage: x = [ k.random_element() for _ in range(500) ]
            sage: y = [ k.random_element() for _ in range(500) ]
            sage: z = KarClass.mul_list(x,y)
            sage: S(z) == S(x)*S(y)
            True
        """
        if self._cutoff < 0:
            return self.mul_step(x,y)
        cdef Py_ssize_t dx = len(x)-1, dy = len(y)-1
        if dx < dy:
            return self.mul_iter(y,x,1)
        else:
            return self.mul_iter(x,y,0)


    cdef list mul_step(self, list x, list y):
        """
        Multiplication step in Karatsuba algorithm
        """
        cdef Py_ssize_t dx = len(x)-1, dy = len(y) - 1
        if dx < 0 or dy < 0: return [ ]
        cdef Py_ssize_t i, j, start, end = dx if dx < self._order else self._order-1
        cdef list twists = [ y ]
        for i from 0 <= i < end:
            twists.append([ self._parent.twist_map()(a) for a in twists[i] ])
        cdef list res = [ ]
        for j from 0 <= j <= dx+dy:
            start = 0 if j <= dy else j-dy
            end = j if j <= dx else dx
            sum = x[start] * twists[start % self._order][j-start]
            for i from start < i <= end:
                sum += x[i] * twists[i % self._order][j-i]
            res.append(sum)
        return res


    cdef list mul_step_matrix(self, list x, list y):
        r"""
        Multiplication step in Karatsuba-matrix algorithm. It
        is based on the ring isomorphism:

        .. MATH::

            S / NS \simeq M_r(k^\sigma)

        with the standard notations::

        -  `S` is the underlying skew polynomial ring

        -  `k` is the base ring of `S` (it's a finite field)

        -  `\sigma` is the twisting automorphism on `k`

        -  `r` is the order of `\sigma`

        -  `N` is a polynomial of definition of the extension
           `k / k^\sigma`

        .. WARNING::

            The polynomials `x` and `y` must have degree
            at most `r^2/2`
        """
        cdef Py_ssize_t dx = len(x)-1, dy = len(y) - 1
        cdef Py_ssize_t i, j
        cdef Py_ssize_t r = self._order
        cdef Mx, My, M
        cdef list row
        cdef RingElement c
        k = self._parent.base_ring()
        Mx = self._T * matrix(k, r, r, x + [self._zero] * (r*r-len(x)))
        My = self._T * matrix(k, r, r, y + [self._zero] * (r*r-len(y)))
        for i from 0 <= i < r:
            frob = self._parent.twist_map(i)
            row = Mx.row(i).list()
            row = map(frob,row)
            row = [ self._t*c for c in row[r-i:] ] + row[:r-i]
            Mx.set_row(i,row)
            row = My.row(i).list()
            row = map(frob,row)
            row = [ self._t*c for c in row[r-i:] ] + row[:r-i]
            My.set_row(i,row)
        M = Mx * My
        for i from 0 <= i < r:
            frob = self._parent.twist_map(-i)
            row = M.row(i).list()
            row = row[i:] + [ c/self._t for c in row[:i] ]
            row = map(frob,row)
            M.set_row(i,row)
        M = self._Tinv * M
        return M.list()[:dx+dy+1]


    cdef list mul_iter(self, list x, list y, char flag): # Assume dx >= dy
        """
        Karatsuba's recursive iteration

        .. WARNING::

            This function assumes that len(x) >= len(y).
        """
        cdef Py_ssize_t i, j, k
        cdef Py_ssize_t dx = len(x)-1, dy = len(y)-1
        if self._algo_matrix:
            if dx < self._cutoff:
                if flag:
                    return self.mul_step_matrix(y,x)
                else:
                    return self.mul_step_matrix(x,y)
        else:
            if dy < self._cutoff:
                if flag:
                    return self.mul_step(y,x)
                else:
                    return self.mul_step(x,y)
        cdef Py_ssize_t dp = self._order * (1 + ceil(dx/self._order/2))
        cdef list x1 = x[:dp], x2 = x[dp:]
        cdef list y1 = y[:dp], y2 = y[dp:]
        cdef list p1, p2, res
        res = self.mul_iter(x1,y1,flag)
        if dy >= dp:
            res.append(self._zero)
            p1 = self.mul_iter(x2,y2,flag)
            res.extend(p1)
            for i from 0 <= i <= dx-dp: x1[i] += x2[i]
            for i from 0 <= i <= dy-dp: y1[i] += y2[i]
            p2 = self.mul_iter(x1,y1,flag)
            j = dx + dy - 2*dp
            k = min(2*dp-2, len(res)-dp-1)
            for i from k >= i > j:
                res[i+dp] += p2[i] - res[i]
            for i from j >= i >= 0:
                res[i+dp] += p2[i] - p1[i] - res[i]
        else:
            if dx-dp < dy:
                p2 = self.mul_iter(y1,x2,not flag)
            else:
                p2 = self.mul_iter(x2,y1,flag)
            for i from 0 <= i < dy:
                res[i+dp] += p2[i]
            res.extend(p2[dy:])
        return res


    def div(self,left,right):
        """
        INPUT:

        -  ``left`` -- a skew polynomial in the skew polynomial
           right attached to the class

        -  ``right`` -- a skew polynomial in the skew polynomial
           right attached to the class

        OUTPUT:

        The quotient and the remainder of the right euclidien division
        of left by right (computed by Karatsuba's algorithm).

        .. WARNING::

            This algorithm is only efficient for very very large
            degrees. Do not use it!

        .. TODO::

            Try to understand why...

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: KarClass = S._karatsuba_class

            sage: a = S.random_element(degree=1000)
            sage: b = S.random_element(degree=500)
            sage: q,r = KarClass.div(a,b)
            sage: q2,r2 = a.quo_rem(b)
            sage: q == q2
            True
            sage: r == r2
            True
        """
        cdef list a = left.list()
        cdef list b = right.list()
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1, i
        cdef list q
        self._twinv = [ ~b[db] ]
        for i from 0 <= i < min(da-db,self._order-1):
            self._twinv.append(self._parent.twist_map()(self._twinv[i]))
        if self._cutoff < 0:
            q = self.div_step(a,0,da,b,0,db)
        else:
            q = self.div_iter(a,0,da,b,0,db)
        return self._parent(q), self._parent(a[:db])


    cdef list div_step(self, list a, Py_ssize_t ia, Py_ssize_t da, list b, Py_ssize_t ib, Py_ssize_t db):
        """
        Division step in Karatsuba's algorithm
        """
        cdef Py_ssize_t i, j
        cdef list q = [ ]
        cdef list twb = [ b[ib:ib+db] ]
        for i from 0 <= i < min(da-db,self._order-1):
            twb.append([ self._parent.twist_map()(x) for x in twb[i] ])
        for i from da-db >= i >= 0:
            c = self._twinv[i%self._order] * a[ia+i+db]
            for j from 0 <= j < db:
                a[ia+i+j] -= c * twb[i%self._order][j]
            q.append(c)
        q.reverse()
        return q


    cdef list div_iter(self, list a, Py_ssize_t ia, Py_ssize_t da, list b, Py_ssize_t ib, Py_ssize_t db):
        """
        Karatsuba recursive iteration
        """
        cdef Py_ssize_t delta = da - db
        if delta < self._cutoff:
            return self.div_step(a,ia,da,b,ib,db)
        cdef Py_ssize_t i
        cdef Py_ssize_t dp = self._order * (1 + int(delta/self._order/2))
        cdef list q = self.div_iter(a,ia+db,delta,b,ib+db-dp,dp), pr
        if db < delta:
            pr = self.mul_iter(q,b[ib:ib+db-dp],0)
        else:
            pr = self.mul_iter(b[ib:ib+db-dp],q,1)
        for i from 0 <= i < len(pr):
            a[i+ia+dp] -= pr[i]
        cdef list qq = self.div_iter(a,ia,db+dp-1,b,ib,db)
        qq.extend(q)
        return qq
