r"""
Univariate Dense Skew Polynomials over a field equipped with a finite order automorphism

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
from sage.functions.other import ceil
from sage.matrix.constructor import zero_matrix, matrix
from sage.rings.ring cimport Ring
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.integer cimport Integer
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.skew_polynomial_element cimport SkewPolynomial_generic_dense

cdef class SkewPolynomial_finite_order_dense(SkewPolynomial_generic_dense):
    def __init__(self, parent, x=None, int check=1, int construct=0, **kwds):
        """
        This method constructs a generic dense skew polynomial over a field equipped
        with an automorphism of finite order.

        INPUT:

        - ``parent`` -- parent of ``self``

        - ``x`` -- list of coefficients from which ``self`` can be constructed

        - ``check`` -- flag variable to normalize the polynomial

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
        SkewPolynomial_generic_dense.__init__ (self, parent, x, check, construct, **kwds)
        self._norm = None
        self._optbound = None


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
        cdef SkewPolynomial_finite_order_dense powx = <SkewPolynomial_finite_order_dense>self._new_c([zero,one], parent)
        cdef SkewPolynomial_finite_order_dense v
        if (exp % 2 == 1):
            v = <SkewPolynomial_finite_order_dense>self._new_c([zero,one], parent)
            if self.degree() == 1:
                v %= self
        else:
            v = <SkewPolynomial_finite_order_dense>self._new_c([one], parent)
        exp = exp >> 1
        n = 1
        while exp != 0:
            powx = (powx.conjugate(n) * powx) % self
            n = n << 1
            if (exp % 2 == 1):
                v = v.conjugate(n)
                v = (v * powx) % self
            exp = exp >> 1
        l = v.list()
        for i from 0 <= i < len(l):
            phir.set_unsafe(i,col,l[i])
        for j from col+1 <= j < d:
            v <<= 1
            v = v.conjugate(1) % self
            for i from 0 <= i <= v.degree():
                phir.set_unsafe(i,j,v._coeffs[i])
        return phir

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


    def reduced_trace(self):
        r"""
        Return the reduced trace of this skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = S.random_element(degree=3,monic=True); a  #random
            x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: T = a.reduced_trace(); T                      #random
            3*(x^3) + 4

        Note that the parent of `T` is the center of the `S`
        (and not `S` itself)::

            sage: T.parent()
            Center of Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5:
            Univariate Polynomial Ring in (x^3) over Finite Field of size 5
            sage: T.parent() == S.center()
            True

        We check that the reduced trace is additive:

            sage: a = S.random_element(degree=5)
            sage: b = S.random_element(degree=7)
            sage: a.reduced_trace() + b.reduced_trace() == (a+b).reduced_trace()
            True
        """
        order = self.parent()._order
        twist_map = self.parent().twist_map()
        coeffs = [ ]
        for i in range(0, self.degree()+1, order):
            tr = c = self._coeffs[i]
            for _ in range(order-1):
                tr = c + twist_map(tr)
            coeffs.append(tr)
        return self.parent().center()(coeffs)


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

            sage: S(N).is_right_divisible_by(a)
            True
            sage: S(N).is_left_divisible_by(a)
            True

        .. NOTE::

            We really need to coerce first `N` into `S`. Otherwise an
            error occurs::

                 sage: N.is_right_divisible_by(a)
                 Traceback (most recent call last):
                 ...
                 AttributeError: 'sage.rings.polynomial.skew_polynomial_element.CenterSkewPolynomial_generic_dense' object has no attribute 'is_right_divisible_by'

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

            sage: S(b).is_right_divisible_by(a)
            True
            sage: S(b).is_left_divisible_by(a)
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

            sage: a = S(Z.random_element(degree=4))
            sage: b = a.bound()
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

            sage: S(b).is_right_divisible_by(a)
            True
            sage: S(b).is_left_divisible_by(a)
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


     # TODO:
     # fast multiplication
     # reduced trace, reduced characteristic polynomial     
