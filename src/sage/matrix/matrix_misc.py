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
                if y.degree() >= d and d >= 0:   x.append(y.coefficients(sparse=False)[d])
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
