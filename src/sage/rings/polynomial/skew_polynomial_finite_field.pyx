r"""
Univariate Dense Skew Polynomials over Finite Fields

This module provides the :class:`~sage.rings.polynomial.skew_polynomial_finite_field.SkewPolynomial_finite_field_dense`
which constructs a single univariate skew polynomial over a finite field equipped with the Frobenius
Endomorphism.

AUTHOR::

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests and refactored classes and methods

"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

import copy
import cysignals
from sage.matrix.constructor import zero_matrix
from sage.rings.ring cimport Ring
from sage.matrix.matrix_dense cimport Matrix_dense
from polynomial_element cimport Polynomial
from sage.rings.integer cimport Integer
from sage.structure.element cimport RingElement
from polynomial_ring_constructor import PolynomialRing
from skew_polynomial_element cimport SkewPolynomial_generic_dense

cdef class SkewPolynomial_finite_field_dense(SkewPolynomial_generic_dense):
    """
    A generic, dense, univariate skew polynomial over a finite field.

    DEFINITION:

    Let `R` be a commutative ring and let `\sigma` be the Frobenius
    Endomorphism over `R` defined by `\sigma(r) = r^p` where `r` is
    an element of `R` and `p` is the prime characteristic of `R`.

    Then, a formal skew polynomial is given by the equation:
    `F(X) = a_{n}X^{n} + ... + a_0`
    where the coefficients `a_{i}` belong to `R` and `X` is a formal variable.

    Addition between two skew polynomials is defined by the usual addition
    operation and the modified multiplication is defined by the rule
    `X a = \sigma(a) X` for all `a` in `R`.

    Skew polynomials are thus non-commutative and the degree of a product
    is equal to the sum of the degree of its factors. The ring of such skew
    polynomials over `R` equipped with `\sigma` is denoted by `S = R[X, \sigma]`
    and it is an additive group.

    .. NOTE::

        #. `S` is a left (resp. right) euclidean noncommutative ring

        #. in particular, every left (resp. right) ideal is principal

    EXAMPLES::

        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]; S
        Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5

    .. SEE ALSO::

        :mod:`sage.rings.polynomial.skew_polynomial_element.SkewPolynomial_generic_dense`
        :mod:`sage.rings.polynomial.skew_polynomial_element.SkewPolynomial`
        :mod:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_field`

    .. TODO::

        Try to replace as possible ``finite field`` by ``field
        endowed with a finite order twist morphism``. It may cause
        new phenomena due to the non trivality of the Brauer group.
    """
    def __init__(self, parent, x=None, int check=1, is_gen=False, int construct=0, **kwds):
        """
        This method constructs a generic dense skew polynomial over a finite field.

        INPUT::

        - ``parent`` -- parent of ``self``

        - ``x`` -- list of coefficients from which ``self`` can be constructed

        - ``check`` -- flag variable to normalize the polynomial

        - ``is_gen`` -- boolean (default: ``False``)

        - ``construct`` -- boolean (default: ``False``)

        TESTS::

            sage: R.<t> = GF(5^3)
            sage: Frob = R.frobenius_endomorphism()
            sage: S.<x> = R['x',Frob]; S
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5

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
        SkewPolynomial_generic_dense.__init__ (self, parent, x, check, is_gen, construct, **kwds)

    cdef SkewPolynomial_finite_field_dense _rgcd(self, SkewPolynomial_finite_field_dense other):
        """
        Fast creation of the right gcd of ``self`` and ``other``.
        """
        cdef SkewPolynomial_finite_field_dense A = self
        cdef SkewPolynomial_finite_field_dense B = other
        cdef SkewPolynomial_finite_field_dense swap
        if len(B._coeffs):
            A = <SkewPolynomial_finite_field_dense>self._new_c(A._coeffs[:],A._parent)
            B = <SkewPolynomial_finite_field_dense>B._new_c(B._coeffs[:],B._parent)
            while len(B._coeffs):
                A._inplace_rrem(B)
                swap = A; A = B; B = swap
            return A
        else:
            return self

    cdef void _inplace_lrem(self, SkewPolynomial_finite_field_dense other):
        """
        Replace ``self`` by the remainder in the left euclidean division
        of ``self`` by ``other`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other)._coeffs
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

    cdef void _inplace_rrem(self, SkewPolynomial_finite_field_dense other):
        """
        Replace ``self`` by the remainder in the right euclidean division
        of ``self`` by ``other`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other)._coeffs
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

    cdef void _inplace_lfloordiv(self, SkewPolynomial_finite_field_dense other):
        """
        Replace ``self`` by the quotient in the left euclidean division
        of ``self`` by ``other`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other)._coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j, deb
        cdef RingElement c, inv
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            (<SkewPolynomial_finite_field_dense>self)._coeffs = [ ]
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

    cdef void _inplace_rfloordiv(self, SkewPolynomial_finite_field_dense other):
        """
        Replace ``self`` by the quotient in the right euclidean division
        of ``self`` by ``other`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other)._coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j, deb, order
        cdef RingElement c, x, inv
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            (<SkewPolynomial_finite_field_dense>self)._coeffs = [ ]
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

    cdef void _inplace_lmonic(self):
        """
        Replace ``self`` by ``self.lmonic()`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef Py_ssize_t da = len(a)-1, i
        cdef RingElement inv = ~a[da]
        parent = self._parent
        a[da] = parent.base_ring()(1)
        for i from 0 <= i < da:
            a[i] *= parent.twist_map(i-da)(inv)

    cdef void _inplace_rmonic(self):
        """
        Replace ``self`` by ``self.rmonic()`` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef Py_ssize_t da = len(a)-1, i
        cdef RingElement inv = ~a[da]
        a[da] = self._parent.base_ring()(1)
        for i from 0 <= i < da:
            a[i] *= inv

    cdef void _inplace_rgcd(self,SkewPolynomial_finite_field_dense other):
        """
        Replace ``self`` by its right gcd with ``other`` (only for internal use).
        """
        cdef SkewPolynomial_finite_field_dense B
        cdef list swap
        if len(other._coeffs):
            B = <SkewPolynomial_finite_field_dense>self._new_c(other._coeffs[:],other._parent)
            while len(B._coeffs):
                B._conjugates = [ B._coeffs ]
                self._inplace_rrem(B)
                swap = self._coeffs
                self._coeffs = B._coeffs
                B._coeffs = swap


    cdef SkewPolynomial_finite_field_dense _rquo_inplace_rem(self, SkewPolynomial_finite_field_dense other):
        """
        Replace ``self`` by the remainder in the right euclidean division
        of ``self`` by ``other`` and return the quotient (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef list b = (<SkewPolynomial_finite_field_dense>other)._coeffs
        cdef Py_ssize_t da = len(a)-1, db = len(b)-1
        cdef Py_ssize_t i, j
        cdef RingElement c, inv
        cdef list q
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError
        if da < db:
            r = self._new_c([],self._parent)
            return r
        inv = ~b[db]
        q = [ ]
        for i from da-db >= i >= 0:
            c = parent.twist_map(i)(inv) * a[i+db]
            q.append(c)
            for j from 0 <= j < db:
                a[i+j] -= c * parent.twist_map(i)(b[j])
        del a[db:]
        self.__normalize()
        q.reverse()
        r = self._new_c(q,self._parent)
        return r

    cdef Py_ssize_t _val_inplace_unit(self):
        """
        Return `v` the valuation of ``self`` and replace ``self`` by
        `self >> v` (only for internal use).
        """
        cdef list a = (<SkewPolynomial_finite_field_dense>self)._coeffs
        cdef Py_ssize_t val = 0
        if len(a) < 0:
            return -1
        while a[0].is_zero():
            del a[0]
            val += 1
        return val

    cdef Matrix_dense _matmul_c(self):
        r"""
        Return the matrix of the multiplication by ``self`` on
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
        if self._norm is None:
            center = self.parent().center()
            if self.is_zero():
                self._norm = center(0)
            else:
                section = center._embed_basering.section()
                exp = (self.parent().base_ring().cardinality() - 1) / (center.base_ring().cardinality() - 1)
                order = self.parent()._order
                lc = section(self.leading_coefficient()**exp)
                if order < self.degree():
                    M = self._matmul_c()
                    self._norm = center([ lc*section(x) for x in M.determinant().monic().list() ])
                else:
                    charpoly = self._matphir_c().characteristic_polynomial()
                    self._norm = center([ lc*section(x) for x in charpoly.list() ])
        return self._norm


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
