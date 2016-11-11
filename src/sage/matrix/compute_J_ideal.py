r"""
Minimal `J`-ideals of matrices

Minimal `J`-ideals of matrices
==============================

Let `B` be an `n\times n`-matrix over a principal ideal domain `D`.

For an ideal `J`, the `J`-ideal of `B` is defined to be
`N_J(B) = \{ f\in D[X] \mid f(B) \in M_n(J) \}`.

For a prime element `p` of `D` and `t\ge 0`, a `(p^t)`-minimal polynomial of `B`
is a monic polynomial `f\in N_{(p^t)}(B)` of minimal degree.

This module computes these minimal polynomials.

Let `p` be a prime element of `D`. Then there is a finite set `\mathcal{S}_p` of
positive integers and monic polynomials `\nu_{ps}` for `s\in\mathcal{S}_p` such
that for `t\ge 1`,

.. MATH::

   N_{(p^t)}(B) = \mu_BD[X] + p^tD[X]
   + \sum_{\substack{s\in\mathcal{S}_p \\ s \le  b(t) }}
   p^{\max\{0,t-s\}}\nu_{ps}D[X]

holds where `b(t) = \min\{r\in \mathcal{S}_p \mid r \ge s\}`. The
degree of `\nu_{ps}` is strictly increasing in `s\in \mathcal{S}_p` and
`\nu_{ps}` is a `(p^s)`-minimal polynomial. If `t\le \max\mathcal{S}_p`,
then the summand `\mu_BD[X]` can be omitted.

TODO: integer valued polynomials

EXAMPLES::

    sage: from calculate_nu import ComputeMinimalPolynomials # not tested
    sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
    sage: C = ComputeMinimalPolynomials(B)
    sage: C.prime_candidates()
    [2, 3, 5]
    sage: for t in range(4):
    ....:     print C.null_ideal(2^t)
    Ideal (1, x^3 + x^2 - 12*x - 20) of Univariate Polynomial
    Ring in x over Integer Ring
    Ideal (2, x^3 + x^2 - 12*x - 20, x^2 + x) of Univariate Polynomial
    Ring in x over Integer Ring
    Ideal (4, x^3 + x^2 - 12*x - 20, x^2 + 3*x + 2) of Univariate Polynomial
    Ring in x over Integer Ring
    Ideal (8, x^3 + x^2 - 12*x - 20, 2*x^2 + 6*x + 4) of Univariate Polynomial
    Ring in x over Integer Ring
    sage: C.p_minimal_polynomials(2)
    {2: x^2 + 3*x + 2}

.. TODO::

   Test code over PIDs other than ZZ.

REFERENCES:

.. [R2016] Roswitha Rissner, Null ideals of matrices over residue class rings of
   principal ideal domains. Linear Algebra Appl., 494:44--69, 2016.

.. [HR2016] Clemens Heuberger and Roswitha Rissner, Computing `J`-Ideals of a
   Matrix Over a Principal Ideal Domain, Preprint, 2016.

AUTHORS:

- Clemens Heuberger (2016)
- Roswitha Rissner (2016)

ACKNOWLEDGEMENT:

- Clemens Heuberger is supported by the Austrian Science Fund (FWF):
  P 24644-N26.

- Roswitha Rissner is supported by the Austrian Science Fund (FWF):
  P 27816-N26.


Classes and Methods
===================
"""

# *****************************************************************************
# Copyright (C) 2016 Clemens Heuberger <clemens.heuberger@aau.at>
#               2016 Roswitha Rissner <roswitha.rissner@tugraz.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

import heapq

from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject
from time import time


def lifting(p, t, A, G):
    r"""
    Compute generators of `\{f \in D[X]^d \mid Af \equiv 0 \pmod{p^{t}}\}` given
    generators of `\{f\in D[X]^d \mid Af \equiv 0\pmod{p^{t-1}}\}`.

    INPUT:

    - ``p`` -- a prime element of some principal ideal domain `D`

    - ``t`` -- a non-negative integer

    - ``A`` -- a `c\times d` matrix  over `D[X]`

    - ``G`` -- a matrix over `D[X]`. The columns of
      `\begin{pmatrix}p^{t-1}I& G\end{pmatrix}` are generators
      of `\{ f\in D[X]^d \mid Af \equiv 0\pmod{p^{t-1}}\}`.
      Can be set to ``None`` if ``t`` is zero.

    OUTPUT:

    A matrix `F` over `D[X]` such that the columns of
    `\begin{pmatrix}p^tI&F&pG\end{pmatrix}` are generators of
    `\{ f\in D[X]^d \mid Af \equiv 0\pmod{p^t}\}`.

    EXAMPLES::

        sage: X = polygen(ZZ, 'X')
        sage: A = matrix([[1, X], [2*X, X^2]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        []
        sage: (A*G1 % 5).is_zero()
        True
        sage: A = matrix([[1, X, X^2], [2*X, X^2, 3*X^3]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        [3*X^2]
        [    X]
        [    1]
        sage: (A*G1 % 5).is_zero()
        True
        sage: G2 = lifting(5, 2, A, G1); G2
        [15*X^2 23*X^2]
        [   5*X      X]
        [     5      1]
        sage: (A*G2 % 25).is_zero()
        True
        sage: lifting(5, 10, A, G1)
        Traceback (most recent call last):
        ...
        ValueError: A*G is not zero mod 5^9

    ALGORITHM:

    [HR2016]_, Algorithm 1.

    TESTS::

        sage: A = matrix([[1, X], [X, X^2]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        Traceback (most recent call last):
        ...
        ValueError: [  1   X|]
        [  X X^2|] does not have full rank
    """
    DX = A.parent().base()
    (X,) = DX.gens()
    D = DX.base_ring()
    d = A.ncols()
    c = A.nrows()

    if t == 0:
        return matrix(DX, d, 0)

    if not (A*G % p**(t-1)).is_zero():
        raise ValueError(
            "A*G is not zero mod %s^%s" % (p, t-1))


    R = A*G/p**(t-1)
    R.change_ring(DX)

    AR = matrix.block([[A, R]])
    Fp = D.quotient(p*D)
    FpX = PolynomialRing(Fp, name=X)

    ARb = AR.change_ring(FpX)
    (Db, Sb, Tb) = ARb.smith_form()
    assert Sb * ARb * Tb == Db
    assert all(i == j or Db[i, j].is_zero()
               for i in range(Db.nrows())
               for j in range(Db.ncols()))

    r = Db.rank()
    if r != c:
        raise ValueError("{} does not have full rank".format(ARb))

    T = Tb.change_ring(DX)

    F1 = matrix.block([[p**(t-1) * matrix.identity(d), G]])*T
    F = F1.matrix_from_columns(range(r, F1.ncols()))
    assert (A*F % (p**t)).is_zero(), "A*F=%s" % (A*F)

    return F


def p_part(f, p):
    r"""
    Compute the `p`-part of a polynomial.

    INPUT:

    - ``f`` -- a polynomial over `D`

    - ``p`` -- a prime in `D`

    OUTPUT:

    A polynomial `g` such that `\deg g \le \deg f` and
    all non-zero coefficients of `f - p g` are not divisible by `p`.

    EXAMPLES::

        sage: X = polygen(ZZ, 'X')
        sage: f = X^3 + 5*X + 25
        sage: g = p_part(f, 5); g
        X + 5
        sage: f - 5*g
        X^3
    """
    DX = f.parent()
    (X,) = DX.gens()
    return sum(c//p * X**i for
               i, c in enumerate(f.list())
               if c % p == 0)


class ComputeMinimalPolynomials(SageObject):
    r"""
    Create an object for computing `(p^t)`-minimal polynomials and `J`-ideals.

    INPUT:

    - ``B`` -- a square matrix over a principal ideal domain `D`.

    OUTPUT:

    An object which allows to call ``p_minimal_polynomials`` and
    ``null_ideal``. TODO: integer valued polynomials

    For an ideal `J`, the `J`-ideal of `B` is defined to be
    `N_J(B) = \{ f\in D[X] \mid f(B) \in M_n(J) \}`.

    For a prime element `p` of `D` and `t\ge 0`, a `(p^t)`-minimal polynomial of
    `B` is a monic polynomial `f\in N_{(p^t)}(B)` of minimal degree.

    The characteristic polynomial of `B` is denoted by `\chi_B`; `n` is the size
    of `B`.

    EXAMPLES::

        sage: from calculate_nu import ComputeMinimalPolynomials # not tested
        sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
        sage: C = ComputeMinimalPolynomials(B)
        sage: for t in range(4):
        ....:     print C.null_ideal(2^t)
        Ideal (1, x^3 + x^2 - 12*x - 20) of Univariate Polynomial
        Ring in x over Integer Ring
        Ideal (2, x^3 + x^2 - 12*x - 20, x^2 + x) of Univariate Polynomial
        Ring in x over Integer Ring
        Ideal (4, x^3 + x^2 - 12*x - 20, x^2 + 3*x + 2) of Univariate Polynomial
        Ring in x over Integer Ring
        Ideal (8, x^3 + x^2 - 12*x - 20, 2*x^2 + 6*x + 4) of Univariate
        Polynomial Ring in x over Integer Ring
        sage: C.p_minimal_polynomials(2)
        {2: x^2 + 3*x + 2}
    """
    def __init__(self, B):
        r"""
        Initialize the ComputeMinimalPolynomials class.

        INPUT:

        - ``B`` -- a square matrix.

        TESTS::

            sage: ComputeMinimalPolynomials(matrix([[1, 2]]))
            Traceback (most recent call last):
            ...
            TypeError: square matrix required.
        """
        super(ComputeMinimalPolynomials, self).__init__()
        if not B.is_square():
            raise TypeError("square matrix required.")

        self._B = B
        X = polygen(B.base_ring())
        adjoint = (X - B).adjoint()
        d = B.nrows()**2
        b = matrix(d, 1, adjoint.list())
        self.chi_B = B.charpoly(X)
        self.mu_B = B.minimal_polynomial()
        self._A = matrix.block([[b , -self.chi_B*matrix.identity(d)]])
        self._DX = X.parent()


    def find_monic_replacements(self, p, t, pt_generators, prev_nu):
        r"""
        Replace possibly non-monic generators of `N_{(p^t)}(B)` by monic
        generators.

        INPUT:

        - ``p`` -- a prime element of `D`

        - ``t`` -- a non-negative integer

        - ``pt_generators`` -- a list `(g_1, \ldots, g_s)` of polynomials in
          `D[X]` such that `N_{(p^t)}(B) = (g_1, \ldots, g_s) + pN_{(p^{t-1})}(B)`.

        - ``prev_nu`` -- a `(p^{t-1})`-minimal polynomial of `B`.

        OUTPUT:

        A list `(h_1, \ldots, h_r)` of monic polynomials such that
        `N_{(p^t)}(B) = (h_1, \ldots, h_r) + pN_{(p^{t-1})}(B)`.

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: x = polygen(ZZ, 'x')
            sage: nu_1 = x^2 + x
            sage: generators_4 = [2*x^2 + 2*x, x^2 + 3*x + 2]
            sage: C.find_monic_replacements(2, 2, generators_4, nu_1)
            [x^2 + 3*x + 2]

        TESTS::

            sage: C.find_monic_replacements(2, 3, generators_4, nu_1)
            Traceback (most recent call last):
            ...
            ValueError: [2*x^2 + 2*x, x^2 + 3*x + 2] are not in N_{(2^3)}(B)
            sage: C.find_monic_replacements(2, 2, generators_4, x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 is not in N_{(2^1)}(B)

        ALGORITHM:

        [HR2016]_, Algorithms 2 and 3.
        """
        if not all((g(self._B) % p**t).is_zero()
                   for g in pt_generators):
            raise ValueError("%s are not in N_{(%s^%s)}(B)" %
                             (pt_generators, p, t))

        if not (prev_nu(self._B) % p**(t-1)).is_zero():
            raise ValueError("%s is not in N_{(%s^%s)}(B)" % (prev_nu, p, t-1))

        (X,) = self._DX.gens()

        replacements = []
        for f in pt_generators:
            g = f
            p_prt = p_part(g, p)

            while g != p*p_prt:
                r = p_prt.quo_rem(prev_nu)[1]
                g1 = g - p*p_prt
                d, u, v = xgcd(g1.leading_coefficient(), p)
                h = u*(p*r + g1) + v*p*prev_nu*X**(g1.degree()-prev_nu.degree())
                replacements.append(h % p**t)
                #reduce coefficients mod p^t to keep coefficients small
                g = g.quo_rem(h)[1]
                p_prt = p_part(g, p)

        replacements = list(set(replacements))
        assert all(g.is_monic() for g in replacements),\
            "Something went wrong in find_monic_replacements"

        return replacements


    def current_nu(self, p, t, pt_generators, prev_nu):
        r"""
        Compute `(p^t)`-minimal polynomial of `B`.

        INPUT:

        - ``p`` -- a prime element of `D`

        - ``t`` -- a positive integer

        - ``pt_generators`` -- a list `(g_1, \ldots, g_s)` of polynomials in
          `D[X]` such that `N_{(p^t)}(B) = (g_1, \ldots, g_s) + pN_{(p^{t-1})}(B)`.

        - ``prev_nu`` -- a `(p^{t-1})`-minimal polynomial of `B`.

        OUTPUT:

        A `(p^t)`-minimal polynomial of `B`.

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: x = polygen(ZZ, 'x')
            sage: nu_1 = x^2 + x
            sage: generators_4 = [2*x^2 + 2*x, x^2 + 3*x + 2]
            sage: C.current_nu(2, 2, generators_4, nu_1)
            x^2 + 3*x + 2

        ALGORITHM:

        [HR2016]_, Algorithm 4.

        TESTS::

            sage: C.current_nu(2, 3, generators_4, nu_1)
            Traceback (most recent call last):
            ...
            ValueError: [2*x^2 + 2*x, x^2 + 3*x + 2] are not in N_{(2^3)}(B)
            sage: C.current_nu(2, 2, generators_4, x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^2 is not in N_{(2^1)}(B)
        """
        if not all((g(self._B) % p**t).is_zero()
                   for g in pt_generators):
            raise ValueError("%s are not in N_{(%s^%s)}(B)" %
                             (pt_generators, p, t))

        if not (prev_nu(self._B) % p**(t-1)).is_zero():
            raise ValueError("%s is not in N_{(%s^%s)}(B)" % (prev_nu, p, t-1))

        generators = self.find_monic_replacements(p, t, pt_generators, prev_nu)

        verbose("------------------------------------------")
        verbose(pt_generators)
        verbose("Generators with (p^t)-generating property:")
        verbose(generators)

        heap = list((f.degree(), f) for f in generators)
        heapq.heapify(heap)

        # find poly of minimal degree
        deg_g, g = heapq.heappop(heap)

        # find nu
        while heap:
            deg_f, f = heapq.heappop(heap)
            #take first element in generators not equal g
            r = (f.quo_rem(g)[1]) % p**t
            if r != 0:
                for h in self.find_monic_replacements(p, t, [r], prev_nu):
                    heapq.heappush(heap, (h.degree(), h))
                if heap and heap[0][0] < deg_g:
                    deg_g, g = heapq.heappushpop(heap, (deg_g, g))

            verbose([g] + [h for (deg_h, h) in heap])

        return g


    def mccoy_column(self, p, t, nu):
        r"""
        Compute matrix for McCoy's criterion.

        INPUT:

        - ``p`` -- a prime element in `D`

        - ``t`` -- a positive integer

        - ``nu`` -- a `(p^t)`-minimal polynomial of `B`

        OUTPUT:

        An `(n^2 + 1) \times 1` matrix `g` with first entry ``nu`` such that
        `\begin{pmatrix}b& -\chi_B I\end{pmatrix}g \equiv 0\pmod{p^t}` where `b`
        consists of the entries of `\operatorname{adj}(X-B)`.

        EXAMPLES::

           sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
           sage: C = ComputeMinimalPolynomials(B)
           sage: x = polygen(ZZ, 'x')
           sage: nu_4 = x^2 + 3*x + 2
           sage: g = C.mccoy_column(2, 2, nu_4)
           sage: b = matrix(9, 1, (x-B).adjoint().list())
           sage: M = matrix.block([[b , -B.charpoly(x)*matrix.identity(9)]])
           sage: (M*g % 4).is_zero()
           True

        ALGORITHM:

        [HR2016]_, Algorithm 5.

        TESTS::

           sage: nu_2 = x^2 + x
           sage: C.mccoy_column(2, 2, nu_2)
           Traceback (most recent call last):
           ...
           ValueError: x^2 + x not in (2^2)-ideal

        """
        if not (nu(self._B) % p**t).is_zero():
            raise ValueError(
                "%s not in (%s^%s)-ideal" % (nu, p, t))

        column = matrix(self._DX, self._A.ncols(), 1,
                        [nu] + [(nu*b).quo_rem(self.chi_B)[0]
                                  for b in self._A[:, 0].list()])

        assert (self._A * column % p**t).is_zero(),\
                                 "McCoy column is not correct"

        return  column


    def p_minimal_polynomials(self, p, s_max=None):
        r"""
        Compute `(p^t)`-minimal polynomials `\nu_t` of `B`.

        INPUT:

        - ``p`` -- a prime in `D`

        - ``s_max`` -- a positive integer (Default: ``None``). If set, only
          `(p^t)`-minimal polynomials for ``t <= s_max`` are computed.

        OUTPUT:

        A dictionary. Keys are a finite subset `\mathcal{S}` of the positive
        integers, the values are the associated `(p^s)`-minimal polynomials
        `\nu_s`, `s\in\mathcal{S}`.

        For `0<t\le \max\mathcal{S}`, a `(p^t)`-minimal polynomial is given by
        `\nu_s` where `s=\min\{ r\in\mathcal{S}\mid r\ge t\}`.  For
        `t>\max\mathcal{S}`, the minimal polynomial of `B` is also a
        `(p^t)`-minimal polynomial.

        If ``s_max`` is set, only those `\nu_s` with ``s <= s_max``
        are returned where ``s_max`` is always included even if it is
        not included in the full set `\mathcal{S}` except when the
        minimal polynomial of `B` is also a ``(p^s_max)``-minimal
        polynomial.

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.p_minimal_polynomials(2)
            {2: x^2 + 3*x + 2}
            sage: set_verbose(1)
            sage: C.p_minimal_polynomials(2)
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            ------------------------------------------
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            p = 2, t = 1:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            Result of lifting:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            F =
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            [x^2 + x]
            [      x]
            [      0]
            [      1]
            [      1]
            [  x + 1]
            [      1]
            [      0]
            [      0]
            [  x + 1]
            verbose 1 (...: calculate_nu.py, current_nu)
            ------------------------------------------
            verbose 1 (...: calculate_nu.py, current_nu)
            (x^2 + x)
            verbose 1 (...: calculate_nu.py, current_nu)
            Generators with (p^t)-generating property:
            verbose 1 (...: calculate_nu.py, current_nu)
            [x^2 + x]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            nu = x^2 + x
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            corresponding columns for G
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            [x^2 + x]
            [  x + 2]
            [      0]
            [      1]
            [      1]
            [  x - 1]
            [     -1]
            [     10]
            [      0]
            [  x + 1]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            ------------------------------------------
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            p = 2, t = 2:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            Result of lifting:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            F =
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            [  2*x^2 + 2*x x^2 + 3*x + 2]
            [          2*x         x + 4]
            [            0             0]
            [            2             1]
            [            2             1]
            [      2*x + 2         x + 1]
            [            2            -1]
            [            0            10]
            [            0             0]
            [      2*x + 2         x + 3]
            verbose 1 (...: calculate_nu.py, current_nu)
            ------------------------------------------
            verbose 1 (...: calculate_nu.py, current_nu)
            (2*x^2 + 2*x, x^2 + 3*x + 2)
            verbose 1 (...: calculate_nu.py, current_nu)
            Generators with (p^t)-generating property:
            verbose 1 (...: calculate_nu.py, current_nu)
            [x^2 + 3*x + 2]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            nu = x^2 + 3*x + 2
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            corresponding columns for G
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            [x^2 + 3*x + 2]
            [        x + 4]
            [            0]
            [            1]
            [            1]
            [        x + 1]
            [           -1]
            [           10]
            [            0]
            [        x + 3]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            ------------------------------------------
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            p = 2, t = 3:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            Result of lifting:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            F =
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            [x^3 + 7*x^2 + 6*x x^3 + 3*x^2 + 2*x]
            [        x^2 + 8*x         x^2 + 4*x]
            [                0                 0]
            [                x             x + 4]
            [            x + 4                 x]
            [    x^2 + 5*x + 4           x^2 + x]
            [           -x + 4                -x]
            [             10*x              10*x]
            [                0                 0]
            [        x^2 + 7*x     x^2 + 3*x + 4]
            verbose 1 (...: calculate_nu.py, current_nu)
            ------------------------------------------
            verbose 1 (...: calculate_nu.py, current_nu)
            (x^3 + 7*x^2 + 6*x, x^3 + 3*x^2 + 2*x)
            verbose 1 (...: calculate_nu.py, current_nu)
            Generators with (p^t)-generating property:
            verbose 1 (...: calculate_nu.py, current_nu)
            [x^3 + 7*x^2 + 6*x, x^3 + 3*x^2 + 2*x]
            verbose 1 (...: calculate_nu.py, current_nu)
            [x^3 + 3*x^2 + 2*x]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials)
            nu = x^3 + 3*x^2 + 2*x
            {2: x^2 + 3*x + 2}
            sage: set_verbose(0)
            sage: C.p_minimal_polynomials(2, s_max=1)
            {1: x^2 + x}
            sage: C.p_minimal_polynomials(2, s_max=2)
            {2: x^2 + 3*x + 2}
            sage: C.p_minimal_polynomials(2, s_max=3)
            {2: x^2 + 3*x + 2}

        ALGORITHM:

        [HR2016]_, Algorithm 5.
        """

        deg_mu = self.mu_B.degree()

        t = 0
        p_min_polys = {}
        nu = self._DX(1)
        d = self._A.ncols()
        G = matrix(self._DX, d, 0)


        while True:
            deg_prev_nu = nu.degree()
            t += 1
            verbose("------------------------------------------")
            verbose("p = %s, t = %s:" % (p, t))

            verbose("Result of lifting:")
            verbose("F =")
            verbose(lifting(p, t, self._A, G))

            nu = self.current_nu(p, t, lifting(p, t, self._A, G)[0], nu)

            verbose("nu = %s" % nu)
            if nu.degree() >= deg_mu:
                return p_min_polys

            if nu.degree() == deg_prev_nu:
                G = G.delete_columns([G.ncols() - 1])
                del p_min_polys[t-1]

            column = self.mccoy_column(p, t, nu) 
            verbose("corresponding columns for G")
            verbose(column)

            G = matrix.block([[p * G, column]])
            p_min_polys[t] = nu

            # allow early stopping for small t
            if t == s_max:
                return p_min_polys


    def null_ideal(self, b=0):
        r"""
        Return the `(b)`-ideal `N_{(b)}(B)=\{f\in D[X] \mid f(B)\in M_n(bD)\}`.

        INPUT:

        - ``b`` -- an element of `D` (Default:  0)

        OUTPUT:

        An ideal in `D[X]`.

        EXAMPLES::

            sage: from calculate_nu import compute_nu # not tested
            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.null_ideal()
            Principal ideal (x^3 + x^2 - 12*x - 20)
            of Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(2)
            Ideal (2, x^3 + x^2 - 12*x - 20, x^2 + x)
            of Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(4)
            Ideal (4, x^3 + x^2 - 12*x - 20, x^2 + 3*x + 2)
            of Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(8)
            Ideal (8, x^3 + x^2 - 12*x - 20, 2*x^2 + 6*x + 4)
            of Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(3)
            Ideal (3, x^3 + x^2 - 12*x - 20)
            of Univariate Polynomial Ring in x over Integer Ring
            sage: C.null_ideal(6)
            Ideal (6, x^3 + x^2 - 12*x - 20, 3*x^2 + 3*x)
            of Univariate Polynomial Ring in x over Integer Ring

        .. TODO::

           Remove minimal polynomial if not required.

        """
        generators = [self.mu_B]

        if b != 0:
            generators = [self._DX(b)] + generators
            for (p, t) in factor(b):
                cofactor = b // p**t
                p_polynomials = self.p_minimal_polynomials(p, t)
                generators += [cofactor*p**(t-s)*nu
                               for s, nu in p_polynomials.iteritems()]

            assert all((g(self._B) % b).is_zero() for g in generators), \
                "Polynomials not in %s-ideal" % (b,)

        return self._DX.ideal(generators)


    def prime_candidates(self):
        r"""
        Determine those primes `p` where `\mu_B` might not be a
        `(p)`-minimal polynomial.

        OUTPUT:

        A list of primes.

        EXAMPLES::

             sage: from calculate_nu import compute_nu # not tested
             sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
             sage: C = ComputeMinimalPolynomials(B)
             sage: C.prime_candidates()
             [2, 3, 5]
             sage: C.p_minimal_polynomials(2)
             {2: x^2 + 3*x + 2}
             sage: C.p_minimal_polynomials(3)
             {}
             sage: C.p_minimal_polynomials(5)
             {}

        This means that `3` and `5` were candidates, but actually, `\mu_B` turns
        out to be a `(3)`-minimal polynomial and a `(5)`-minimal polynomial.
        """
        F, T = self._B.frobenius(2)

        return [p for (p, t) in factor(T.det())]
