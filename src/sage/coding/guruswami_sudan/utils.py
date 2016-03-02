r"""
Guruswami-Sudan utility methods

AUTHORS:

- Johan S. R. Nielsen, original implementation (see [Nielsen]_ for details)
- David Lucas, ported the original implementation in Sage
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#                     2015 Johan S. R. Nielsen <jsrn@jsrn.dk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.functions.other import binomial, floor, sqrt
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.arith.all import lcm
from sage.combinat.permutation import Permutation

def polynomial_to_list(p, len):
    r"""
    Returns ``p`` as a list of its coefficients of length ``len``.

    INPUT:

    - ``p`` -- a polynomial

    - ``len`` -- an integer. If ``len`` is smaller than the degree of ``p``, the
      returned list will be of size degree of ``p``, else it will be of size ``len``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import polynomial_to_list
        sage: F.<x> = GF(41)[]
        sage: p = 9*x^2 + 8*x + 37
        sage: polynomial_to_list(p, 4)
        [37, 8, 9, 0]
    """
    return list(p) + [0]*max(0, len-p.degree()-1)

def johnson_radius(n, d):
    r"""
    Returns the Johnson-radius for the code length `n` and the minimum distance `d`.

    The Johnson radius is defined as `n - \sqrt(n(n-d))`.

    INPUT:

    - ``n`` -- an integer, the length of the code
    - ``d`` -- an integer, the minimum distance of the code

    EXAMPLES::

        sage: sage.coding.guruswami_sudan.utils.johnson_radius(250, 181)
        -5*sqrt(690) + 250
    """
    return n - sqrt(n*(n-d))

def ligt(x):
    r"""
    Returns the least integer greater than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import ligt
        sage: ligt(41)
        42

    It works with any type of numbers (not only integers)::

        sage: ligt(41.041)
        42
    """
    return floor(x+1)

def gilt(x):
    r"""
    Returns the greatest integer smaller than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import gilt
        sage: gilt(43)
        42

    It works with any type of numbers (not only integers)::

        sage: gilt(43.041)
        43
    """
    if x in ZZ:
        return Integer(x-1)
    else:
        return floor(x)

def solve_degree2_to_integer_range(a,b,c):
    r"""
    Returns the greatest integer range `[i_1, i_2]` such that
    `i_1 > x_1` and `i_2 < x_2` where `x_1, x_2` are the two zeroes of the equation in `x`:
    `ax^2+bx+c=0`.

    If there is no real solution to the equation, it returns an empty range with negative coefficients.

    INPUT:

    - ``a``, ``b`` and ``c`` -- coefficients of a second degree equation, ``a`` being the coefficient of
      the higher degree term.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import solve_degree2_to_integer_range
        sage: solve_degree2_to_integer_range(1, -5, 1)
        (1, 4)

    If there is no real solution::

        sage: solve_degree2_to_integer_range(50, 5, 42)
        (-2, -1)
    """
    D = b**2 - 4*a*c
    if D < 0:
        return (-2,-1)
    sD = float(sqrt(D))
    minx, maxx = (-b-sD)/2.0/a , (-b+sD)/2.0/a
    mini, maxi = (ligt(minx), gilt(maxx))
    if mini > maxi:
        return (-2,-1)
    else:
        return (mini,maxi)

def apply_weights(M, weights):
    r"""
    Applies column weights inplace to the matrix `M`.

    If ``weights`` are all integers, then `M` is multiplied on the `n`th
    column with `x^{weights[n]}`.

    If weights are fractions, then `M` is appropriately column permuted and
    multiplied with the `x^t` where `t = int(weights[n])`. Afterwards, the
    permutation is returned if needed for reference.

    INPUT:

    - ``M`` -- a matrix

    - ``weights`` -- a list

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import apply_weights
        sage: F.<x> = GF(7)[]
        sage: M = matrix(F, [[2*x^2 + x, 5*x^2 + 2*x + 1, 4*x^2 + x],\
                             [x^2 + 3*x + 3, 5*x^2 + 5*x + 1, 6*x^2 + 5*x + 4],\
                             [5*x^2 + 2*x + 4, 4*x^2 + 2*x, 5*x^2 + x + 2]])
        sage: weights = [1, 2, 3]
        sage: apply_weights(M, weights)
        sage: M
        [          2*x^3 + x^2   5*x^4 + 2*x^3 + x^2           4*x^5 + x^4]
        [    x^3 + 3*x^2 + 3*x   5*x^4 + 5*x^3 + x^2 6*x^5 + 5*x^4 + 4*x^3]
        [  5*x^3 + 2*x^2 + 4*x         4*x^4 + 2*x^3   5*x^5 + x^4 + 2*x^3]
    """
    x = M.base_ring().gen()
    if all(w.is_integer() for w in weights):
        for j in range(M.ncols()):
            M.set_col_to_multiple_of_col(j,j, x**weights[j])
    else:
        perm = fractional_weight_permutation(weights)
        for j in range(M.ncols()):
            M.set_col_to_multiple_of_col(j,j, x**floor(weights[j]))
        M.permute_columns(perm)
        return perm

def remove_weights(M, weights):
    r"""
    Removes the weights inplace to the matrix ``M``
    as they were introduced by :func:`apply_weights`.

    INPUT:

    - ``M`` -- a matrix

    - ``weights`` -- a list

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import remove_weights
        sage: F.<x> = GF(7)[]
        sage: M = matrix(F, [[2*x^3 + x^2, 5*x^4 + 2*x^3 + x^2, 4*x^5 + x^4],\
                            [x^3 + 3*x^2 + 3*x, 5*x^4 + 5*x^3 + x^2, 6*x^5 + 5*x^4 + 4*x^3],\
                            [5*x^3 + 2*x^2 + 4*x, 4*x^4 + 2*x^3, 5*x^5 + x^4 + 2*x^3]])

        sage: weights = [1, 2, 3]
        sage: remove_weights(M, weights)
        sage: M
        [      2*x^2 + x 5*x^2 + 2*x + 1       4*x^2 + x]
        [  x^2 + 3*x + 3 5*x^2 + 5*x + 1 6*x^2 + 5*x + 4]
        [5*x^2 + 2*x + 4     4*x^2 + 2*x   5*x^2 + x + 2]
    """
    if all(w.is_integer() for w in weights):
        for i in range(M.nrows()):
            for j in range(M.ncols()):
                M[i,j] = M[i,j].shift(-weights[j])
    else:
        perm = fractional_weight_permutation(weights)
        pinv = perm.inverse()
        M.permute_columns(pinv)
        remove_weights(M, [floor(wj) for wj in weights])

def fractional_weight_permutation(weights):
    r"""
    Returns the permutation which can be used for embedding the semantics of
    fractional weights into the module minimisation.
    A permutation is returned of the integers from 1 to ``len(numerators)``.

    INPUT:

    - ``weights`` -- a list of fractions

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import fractional_weight_permutation
        sage: weights = [1/4, 1/2, 3/4]
        sage: fractional_weight_permutation(weights)
        [1, 2, 3]
    """
    from sage.coding.guruswami_sudan.interpolation import _flatten_once
    n = len(weights)
    denominator = lcm(list(f.denominator() for f in weights))
    numerators = [f.numerator() * denominator/f.denominator() for f in weights]
    residues = [num % denominator for num in numerators]
    res_map = dict()
    for i in range(n):
        if residues[i] in res_map:
            res_map[residues[i]].append(i+1)
        else:
            res_map[residues[i]] = [i+1]
    res_uniq = sorted(res_map.keys())
    return Permutation(list(_flatten_once([ res_map[res] for res in res_uniq])))

def _leading_position(v, weights=None):
    r"""
    Returns the position of the highest-degree term of ``v``.

    This methods can manage weighted degree, by providing ``weight`` to it.

    In case of several positions having the same, highest degree, the highest position is given.

    INPUT:

    - ``v`` -- a vector of polynomials

    - ``weights`` -- (default: ``None``) a vector of integers or fractions, the weights of ``v``.
      If ``None``, all weights are considered as ``0``.

    EXAMPLES::

    sage: from sage.coding.guruswami_sudan.utils import _leading_position
    sage: F.<x> = GF(7)[]
    sage: v = vector(F, [3*x^2 + 3*x + 4, 4*x + 3, 4*x^2 + 4*x + 5, x^2 + 2*x + 5, 3*x^2 + 5*x])
    sage: _leading_position(v)
    4
    """
    if not weights:
        weights=[0]*len(v)
    best=-1
    bestp=-1
    for p in range(0,len(v)):
        if not v[p].is_zero():
            vpdeg = v[p].degree() + weights[p]
            if vpdeg >= best:
                best=vpdeg
                bestp = p
    if best==-1:
        return -1
    else:
        return bestp

def leading_term(v, weights=None):
    r"""
    Returns the term of ``v`` with the highest degree.

    This methods can manage weighted degree, by providing ``weight`` to it.

    In case of several positions having the same, highest degree, the term with
    the highest position is returned.

    INPUT:

    - ``v`` -- a vector of polynomials

    - ``weights`` -- (default: ``None``) a vector of integers or fractions, the weights of ``v``.
      If ``None``, all weights are considered as ``0``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import leading_term
        sage: F.<x> = GF(7)[]
        sage: v = vector(F, [3*x^2 + 3*x + 4, 4*x + 3, 4*x^2 + 4*x + 5, x^2 + 2*x + 5, 3*x^2 + 5*x])
        sage: leading_term(v)
        3*x^2 + 5*x
    """
    return v[_leading_position(v, weights=weights)]
