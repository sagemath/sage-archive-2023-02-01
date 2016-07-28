r"""
Univariate Dense Skew Polynomials over Finite Fields

This module provides the :class:`~sage.rings.polynomial.skew_polynomial_finite_field.SkewPolynomial_finite_field_dense`
which constructs a single univariate skew polynomial over a finite field equipped with the Frobenius
Endomorphism.

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

import copy
import cysignals
from sage.matrix.constructor import zero_matrix
from sage.rings.ring cimport Ring
from polynomial_ring_constructor import PolynomialRing
from skew_polynomial_element cimport SkewPolynomial_generic_dense

cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_generic_dense):
    """
    A generic, dense, univariate skew polynomial over a finite field.

    DEFINITION:

    Let `R` be a commutative ring and let `\sigma` be the Frobenius
    Endomorphism over `R` defined by:
    `\sigma(r) = r^p`
    where `r` is an element of `R` and `p` is the prime characteristic of `R`.

    Then, a formal skew polynomial is given by the equation:
    `F(X) = a_{n}X^{n} + ... + a_0`
    where the coefficients `a_{i}` belong to `R` and `X` is a formal variable.

    Addition between two skew polynomials is defined by the usual addition
    operation and the modified multiplication is defined by the rule
    `X a = \sigma(a) X` for all `a` in `R`. Skew polynomials are thus
    non-commutative and the degree of a product is equal to the sum of the
    degrees of the factors.

    The ring of such skew polynomials over `R` equipped with `\sigma` is
    denoted by `S = k[X,\sigma]`and it is an addtive group.

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
            sig_off()
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
            sig_off()
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
