"""
Miscellaneous matrix functions

"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.fields import Fields
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
_Fields = Fields()

def row_iterator(A):
    for i in xrange(A.nrows()):
        yield A.row(i)

def weak_popov_form(M,ascend=True):
    """
    This function computes a weak Popov form of a matrix over a rational
    function field `k(x)`, for `k` a field.

    INPUT:

     - `M` - matrix

     - `ascend` - if True, rows of output matrix `W` are sorted so
       degree (= the maximum of the degrees of the elements in
       the row) increases monotonically, and otherwise degrees decrease.

    OUTPUT:

    A 3-tuple `(W,N,d)` consisting of two matrices over `k(x)` and a list
    of integers:

    1. `W` - matrix giving a weak the Popov form of M
    2. `N` - matrix representing row operations used to transform
       `M` to `W`
    3. `d` - degree of respective columns of W; the degree of a column is
       the maximum of the degree of its elements

    `N` is invertible over `k(x)`. These matrices satisfy the relation
    `N*M = W`.

    EXAMPLES:

    The routine expects matrices over the rational function field, but
    other examples below show how one can provide matrices over the ring
    of polynomials (whose quotient field is the rational function field).

    ::

        sage: R.<t> = GF(3)['t']
        sage: K = FractionField(R)
        sage: import sage.matrix.matrix_misc
        sage: sage.matrix.matrix_misc.weak_popov_form(matrix([[(t-1)^2/t],[(t-1)]]))
        (
        [          0]  [      t 2*t + 1]
        [(2*t + 1)/t], [      1       2], [-Infinity, 0]
        )

    NOTES:

    See docstring for weak_popov_form method of matrices for
    more information.
    """

    # determine whether M has polynomial or rational function coefficients
    R0 = M.base_ring()

    #Compute the base polynomial ring
    if R0 in _Fields:
        R = R0.base()
    else:
        R = R0
    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
    if not is_PolynomialRing(R):
        raise TypeError("the coefficients of M must lie in a univariate polynomial ring")

    t = R.gen()

    # calculate least-common denominator of matrix entries and clear
    # denominators. The result lies in R
    from sage.rings.arith import lcm
    from sage.matrix.constructor import matrix
    from sage.misc.functional import numerator
    if R0 in _Fields:
        den = lcm([a.denominator() for a in M.list()])
        num = matrix([(lambda x : map(numerator,  x))(v) for v in map(list,(M*den).rows())])
    else:
        # No need to clear denominators
        den = R.one_element()
        num = M

    r = [list(v) for v in num.rows()]

    N = matrix(num.nrows(), num.nrows(), R(1)).rows()

    from sage.rings.infinity import Infinity
    if M.is_zero():
        return (M, matrix(N), [-Infinity for i in range(num.nrows())])

    rank = 0
    num_zero = 0
    while rank != len(r) - num_zero:
        # construct matrix of leading coefficients
        v = []
        for w in map(list, r):
            # calculate degree of row (= max of degree of entries)
            d = max([e.numerator().degree() for e in w])

            # extract leading coefficients from current row
            x = []
            for y in w:
                if y.degree() >= d and d >= 0:   x.append(y.coeffs()[d])
                else:                            x.append(0)
            v.append(x)
        l = matrix(v)

        # count number of zero rows in leading coefficient matrix
        # because they do *not* contribute interesting relations
        num_zero = 0
        for v in l.rows():
            is_zero = 1
            for w in v:
                if w != 0:
                    is_zero = 0
            if is_zero == 1:
                num_zero += 1

        # find non-trivial relations among the columns of the
        # leading coefficient matrix
        kern = l.kernel().basis()
        rank = num.nrows() - len(kern)

        # do a row operation if there's a non-trivial relation
        if not rank == len(r) - num_zero:
            for rel in kern:
                # find the row of num involved in the relation and of
                # maximal degree
                indices = []
                degrees = []
                for i in range(len(rel)):
                    if rel[i] != 0:
                        indices.append(i)
                        degrees.append(max([e.degree() for e in r[i]]))

                # find maximum degree among rows involved in relation
                max_deg = max(degrees)

                # check if relation involves non-zero rows
                if max_deg != -1:
                    i = degrees.index(max_deg)
                    rel /= rel[indices[i]]

                    for j in range(len(indices)):
                        if j != i:
                            # do row operation and record it
                            v = []
                            for k in range(len(r[indices[i]])):
                                v.append(r[indices[i]][k] + rel[indices[j]] * t**(max_deg-degrees[j]) * r[indices[j]][k])
                            r[indices[i]] = v

                            v = []
                            for k in range(len(N[indices[i]])):
                                v.append(N[indices[i]][k] + rel[indices[j]] * t**(max_deg-degrees[j]) * N[indices[j]][k])
                            N[indices[i]] = v

                    # remaining relations (if any) are no longer valid,
                    # so continue onto next step of algorithm
                    break

    # sort the rows in order of degree
    d = []
    from sage.rings.all import infinity
    for i in range(len(r)):
        d.append(max([e.degree() for e in r[i]]))
        if d[i] < 0:
            d[i] = -infinity
        else:
            d[i] -= den.degree()

    for i in range(len(r)):
        for j in range(i+1,len(r)):
            if (ascend and d[i] > d[j]) or (not ascend and d[i] < d[j]):
                (r[i], r[j]) = (r[j], r[i])
                (d[i], d[j]) = (d[j], d[i])
                (N[i], N[j]) = (N[j], N[i])

    # return reduced matrix and operations matrix
    return (matrix(r)/den, matrix(N), d)

def prm_mul(p1, p2, free_vars_indices):
    """
    Return the product of ``p1`` and ``p2``, putting free variables in
    ``free_vars_indices`` to `1`.

    This function is mainly use as a subroutine of
    :func:`permanental_minor_vector`.

    EXAMPLES::

        sage: from sage.matrix.matrix_misc import prm_mul
        sage: t = polygen(ZZ, 't')
        sage: p1 = {0: 1, 1: t, 4: t}
        sage: p2 = {0: 1, 1: t, 2: t}
        sage: prm_mul(p1, p2, [0])
        {0: 2*t + 1, 2: t^2 + t, 4: t^2 + t, 6: t^2}
    """
    p = {}
    mask_free = 0
    one = 1
    for i in free_vars_indices:
        mask_free += one << i
    if not p2:
        return p
    get = p.get
    for exp1, v1 in p1.iteritems():
        for exp2, v2 in p2.iteritems():
            if exp1 & exp2:
                continue
            exp = exp1 | exp2
            v = v1 * v2
            if exp & mask_free:
                for i in free_vars_indices:
                    if exp & (one << i):
                        exp ^= one << i
            if exp not in p:
                p[exp] = v
            else:
                p[exp] += v
    return p

def permanental_minor_vector(m, permanent_only=False):
    r"""
    Return the polynomial of the sums of permanental minors of a matrix `m`
    in array form.

    The polynomial of the sums of permanental minors is

    .. MATH::

        \sum_0^{min(nrows, ncols)} p_i(m) x^i

    where `p_i(m)` is the `i`-th permanental minor of `m` (that can also be
    obtained through the method
    :meth:`~sage.matrix.matrix2.Matrix.permanental_minor` via
    ``m.permanental_minor(i)``).

    INPUT:

     - ``m`` -- matrix

     - ``permanent_only`` -- optional boolean. If ``True``, only the permanent
       is computed (might be faster).

    OUTPUT:

    The list of coefficients of the polynomial, i.e. the coefficient of `x^i` is
    at the `i`-th position in the list. In particular, the last element of the
    list is the permanent.

    EXAMPLES::

        sage: from sage.matrix.matrix_misc import permanental_minor_vector
        sage: m = matrix([[1,1],[1,2]])
        sage: permanental_minor_vector(m)
        [1, 5, 3]
        sage: permanental_minor_vector(m, permanent_only=1)
        3

    ::

        sage: M = MatrixSpace(ZZ,4,4)
        sage: A = M([1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
        sage: permanental_minor_vector(A)
        [1, 28, 114, 84, 0]
        sage: [A.permanental_minor(i) for i in range(5)]
        [1, 28, 114, 84, 0]

    An example over `\QQ`::

        sage: M = MatrixSpace(QQ,2,2)
        sage: A = M([1/5,2/7,3/2,4/5])
        sage: permanental_minor_vector(A, True)
        103/175

    An example with polynomial coefficients::

        sage: R.<a> = PolynomialRing(ZZ)
        sage: A = MatrixSpace(R,2)([[a,1], [a,a+1]])
        sage: permanental_minor_vector(A, True)
        a^2 + 2*a

    ALGORITHM:

        The permanent of a n x n matrix `A` is the coefficient of the
        x_1...x_n monomial in

        .. MATH::

            \prod_{i=1}^n \sum_{j=1}^n A_{ij} x_j

        Evaluating this product one can neglect `x_i^2`, that is `x_i`
        can be considered to be nilpotent of order `2`.

        To formalize this procedure, consider polynomials over the algebra
        `K[t, \eta_1, \eta_2,\ldots, \eta_k]` where the `\eta_i` are
        commuting, nilpotent of order `2` (i.e. `\eta_i^2 = 0`),
        and `t` is commuting.
        Introduce an "integration" operation <p> consisting in the sum
        of the coefficients of the non-vanishing monomials of the polynomial p.

        Let us consider an example of computation.
        Let `p_1 = 1 + t \eta_0 + t \eta_1` and
        `p_2 = 1 + t \eta_0 + t \eta_2`. Then

        .. MATH::

            p_1 p_2 = 1 + 2t \eta_0 +
                    t (\eta_1 + \eta_2) +
                    t^2 (\eta_0 \eta_1 + \eta_0 \eta_2 + \eta_1 \eta_2)

            <p_1 p_2> = 1 + 4t + 3t^2

        A useful property is the following: let `p_1(\eta_0,..,\eta_k)`
        be a polynomial in distributed form, and let `p_2(\eta_j,..,\eta_k)`
        be a product of polynomials, with `j \ge 1`; one has

        .. MATH::
            <p_1(\eta_0,..,\eta_k) p_2> = <p_1(1,..,\eta_j,..,\eta_k) p_2>

        where `\eta_0,..,\eta_{j-1}` are replaced by `1` in `p_1`.

        In this formalism

        .. MATH::

            perm(A) = < \prod_{i=1}^n \sum_{j=1}^n A_{ij} \eta_j >

        The sum of permanental k-minors of `A` is

        .. MATH::

            permMinor (A, k) = \sum_{r,s} perm(A_{r,s})

        where ``A_{r,s}`` is the matrix obtained eliminating the rows `r`
        and the columns `s`.

        The generating function of the sums of the permanental minors
        of a m x n matrix is

        .. MATH::

            g(t) = < \prod_{i=1}^m (1 + t \sum_{j=1}^n A_{ij} \eta_j) >

        In fact the `t^k` coefficient of `g(t)` corresponds to choosing
        `k` rows of `A`;  `\eta_i` is associated to the i-th column;
        nilpotency avoids having twice the same column in a product of A's.

        The product is implemented as a subroutine in :func:`prm_mul`. The
        polynomials are represented in dictionary form: to a variable `\eta_i`
        it is associated the key `2^i` (or in Python ``1 << i``).  So in the
        above example `\eta_1` corresponds to the key `2` while the product
        `\eta_1 \eta_2` to the key `6`.

        The complexity of this algorithm is `O(2^n m^2 n)`.
        If `A` is a banded matrix with width `w`, that is `A_{ij}=0` for
        `|i - j| > w`, and `w < n/2`, then the complexity of the
        algorithm is `O(4^w (w+1) n^2)`.

    REFERENCES:

    .. [ButPer] P. Butera and M. Pernici "Sums of permanental minors
       using Grassmann algebra", :arxiv:`1406.5337`
    """
    K = PolynomialRing(m.base_ring(), 't')
    nrows = m.nrows()
    ncols = m.ncols()
    m = m.rows()
    p = {0: K.one()}
    t = K.gen()
    done_vars = set()
    one = 1
    for i in range(nrows):
        if permanent_only:
            p1 = {}
        else:
            p1 = {0: K.one()}
        a = m[i]
        for j in range(len(a)):
            if a[j]:
                p1[one<<j] = a[j] * t
        free_vars_indices = []
        for j in range(ncols):
            if j in done_vars:
                continue
            if all(m[k][j] == 0 for k in range(i+1, len(m))):
                free_vars_indices.append(j)
                done_vars.add(j)
        p = prm_mul(p, p1, free_vars_indices)
    if not p:
        return K.zero()
    assert len(p) == 1
    a = p[0].coeffs()
    len_a = min(nrows, ncols) + 1
    if len(a) < len_a:
        a += [0] * (len_a - len(a))
    if permanent_only:
        return a[-1]
    return a
