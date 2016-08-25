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
from sage.functions.other import ceil
from sage.matrix.constructor import zero_matrix, matrix
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
