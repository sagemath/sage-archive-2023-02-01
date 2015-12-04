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

IMPOSSIBLE_PARAMS = "Impossible parameters for the Guruswami-Sudan algorithm"

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

def list_decoding_range(n, d, q=None):
    r"""
    Returns the minimal and maximal number of errors correctable by a
    Johnson-distance list decoder beyond half the minimal distance.

    INPUT:

    - ``n`` -- an integer, the length of the code
    - ``d`` -- an integer, the minimum distance of the code
    - ``q`` -- (default: ``None``) an integer, the field characteristic

    EXAMPLES::

        sage: sage.coding.guruswami_sudan.utils.list_decoding_range(250, 181)
        (91, 118)
    """
    if q is None:
        return (ligt((d-1)/2), gilt(n - sqrt(n*(n-d))))
    else:
        f = (q-1.)/q
        return (ligt((d-1)/2), gilt(f*(n-sqrt(n*(n-d/f)))))

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

    If there is no real solution to the equation, it returns an empty, range with negative coefficients.

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

def find_minimal_satisfiable(f, startn=1, contiguous=True):
    r"""
    Returns the minimal integral ``n``, `n > 0` such that ``f(n) == True``.

    `startn` can be given as a hint to a value that might be true.

    INPUT:

    - ``f`` -- a function

    - ``startn`` -- (default: ``1``) the starting point of the algorithm.

    EXAMPLES::

        sage: def f(x):
        ....:    return None if x > 10 or x == 1 else x + 1

        sage: sage.coding.guruswami_sudan.utils.find_minimal_satisfiable(f)
        2
    """
    maxn = startn
    minn = 1
    # Keep doubling n to find one that works and then switch to binary search
    while not f(maxn):
        minn = maxn + 1
        maxn *= 2
    while minn < maxn:
        tryn = minn + floor((maxn - minn) * 0.5)
        if f(tryn):
            maxn = tryn
        else:
            minn = tryn + 1
    return maxn

def module_apply_weights(M, weights):
    """Apply column weights inplace to the matrix `M`. If weights are all integers, then
    `M` is multiplied on the $n$th column with `x^{weights[n]}`.
    If weights are fractions, then `M` is appropriately column permuted and
    multiplied with the `x^t` where `t = int(weights[n])`. Afterwards, the
    permutation is returned if needed for reference. See e.g. [Dec]
    """
    if all( w.is_integer() for w in weights):
        for j in range(M.ncols()):
            M.set_col_to_multiple_of_col(j,j, x^weights[j])
    else:
        perm = module_fractional_weight_permutation(weights)
        for j in range(M.ncols()):
            M.set_col_to_multiple_of_col(j,j, x^floor(weights[j]))
        M.permute_columns(perm)
        return perm

def module_remove_weights(M, weights):
    """Remove the weights as introduced by ``module_apply_weights``."""
    #M = M.change_ring(QQ)
    if all( w.is_integer() for w in weights):
        for i in range(M.nrows()):
            for j in range(M.ncols()):
                M[i,j] = M[i,j].shift(-weights[j])
    else:
        perm = module_fractional_weight_permutation(weights)
        pinv = perm.inverse()
        M.permute_columns(pinv)
        module_remove_weights(M, [ floor(wj) for wj in weights ])

def module_fractional_weight_permutation(weights):
    """Return the permutation which can be used for embedding the semantics of
    fractional weights into the module minimisation.
    A permutation is returned of the integers from 1 to ``len(numerators)``.
    See e.g. [Dec].
    """
    n = len(weights)
    denominator = lcm(list(f.denominator() for f in weights))
    numerators = [ f.numerator() * denominator/f.denominator() for f in weights ]
    residues = [ num % denominator for num in numerators ]
    res_map = dict()
    for i in range(n):
        if residues[i] in res_map:
            res_map[residues[i]].append(i+1)
        else:
            res_map[residues[i]] = [i+1]
    res_uniq = sorted(res_map.keys())
    return Permutation(list(flatten_once([ res_map[res] for res in res_uniq ])))

