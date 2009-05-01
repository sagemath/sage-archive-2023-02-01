"""
The Victor Miller Basis
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import math

from sage.matrix.all import MatrixSpace, Matrix
from sage.rings.all import QQ, ZZ, Integer, binomial, PowerSeriesRing, O as bigO
from sage.structure.all import Sequence
from sage.libs.all import ntl
from sage.misc.all import verbose

from eis_series import eisenstein_series_qexp

def victor_miller_basis(k, prec=10, cusp_only=False, var='q'):
    r"""
    Compute and return the Victor-Miller basis for modular forms of
    weight k and level 1 to precision $O(q^prec)$.  if
    ``cusp_only`` is True, return only a basis for the cuspidal
    subspace.

    INPUT:

        - ``k`` -- an integer

        - ``prec`` -- (default: 10) a positive integer

        - ``cusp_only`` -- bool (default: False)

        - ``var`` -- string (default: 'q'

    OUTPUT:
        A sequence whose entries are power series in ZZ[[var]].

    EXAMPLES::

        sage: victor_miller_basis(1, 6)
        []
        sage: victor_miller_basis(0, 6)
        [
        1 + O(q^6)
        ]
        sage: victor_miller_basis(2, 6)
        []
        sage: victor_miller_basis(4, 6)
        [
        1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
        ]

        sage: victor_miller_basis(6, 6, var='w')
        [
        1 - 504*w - 16632*w^2 - 122976*w^3 - 532728*w^4 - 1575504*w^5 + O(w^6)
        ]

        sage: victor_miller_basis(6, 6)
        [
        1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
        ]
        sage: victor_miller_basis(12, 6)
        [
        1 + 196560*q^2 + 16773120*q^3 + 398034000*q^4 + 4629381120*q^5 + O(q^6),
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
        ]

        sage: victor_miller_basis(12, 6, cusp_only=True)
        [
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
        ]
        sage: victor_miller_basis(24, 6, cusp_only=True)
        [
        q + 195660*q^3 + 12080128*q^4 + 44656110*q^5 + O(q^6),
        q^2 - 48*q^3 + 1080*q^4 - 15040*q^5 + O(q^6)
        ]
        sage: victor_miller_basis(24, 6)
        [
        1 + 52416000*q^3 + 39007332000*q^4 + 6609020221440*q^5 + O(q^6),
        q + 195660*q^3 + 12080128*q^4 + 44656110*q^5 + O(q^6),
        q^2 - 48*q^3 + 1080*q^4 - 15040*q^5 + O(q^6)
        ]
        sage: victor_miller_basis(32, 6)
        [
        1 + 2611200*q^3 + 19524758400*q^4 + 19715347537920*q^5 + O(q^6),
        q + 50220*q^3 + 87866368*q^4 + 18647219790*q^5 + O(q^6),
        q^2 + 432*q^3 + 39960*q^4 - 1418560*q^5 + O(q^6)
        ]

        sage: victor_miller_basis(40,200)[1:] == victor_miller_basis(40,200,cusp_only=True)
        True
        sage: victor_miller_basis(200,40)[1:] == victor_miller_basis(200,40,cusp_only=True)
        True
    """
    k = Integer(k)
    if k%2 == 1 or k==2:
        return Sequence([])
    elif k < 0:
        raise ValueError, "k must be non-negative"
    elif k == 0:
        return Sequence([PowerSeriesRing(ZZ,var)(1).add_bigoh(prec)], cr=True)
    e = k.mod(12)
    if e == 2: e += 12
    n = (k-e) // 12

    if n == 0 and cusp_only:
        return Sequence([])

    # If prec is less than or equal to the dimension of the space of
    # cusp forms, which is just n, then we know the answer, and we
    # simply return it.
    if prec <= n:
        q = PowerSeriesRing(ZZ,var).gen(0)
        err = bigO(q**prec)
        ls = [0] * (n+1)
        if not cusp_only:
            ls[0] = 1 + err
        for i in range(1,prec):
            ls[i] = q**i + err
        for i in range(prec,n+1):
            ls[i] = err
        return Sequence(ls, cr=True)

    F6 = ((-504)*eisenstein_series_qexp(6,prec,var=var)).change_ring(ZZ)

    if e == 0:
        A = F6.parent()(1)
    elif e == 4:
        A = (240*eisenstein_series_qexp(4,prec,var=var)).change_ring(ZZ)
    elif e == 6:
        A = F6
    elif e == 8:
        A = (480*eisenstein_series_qexp(8,prec,var=var)).change_ring(ZZ)
    elif e == 10:
        A = (-264*eisenstein_series_qexp(10,prec,var=var)).change_ring(ZZ)
    else: # e == 14
        A = (-24*eisenstein_series_qexp(14,prec,var=var)).change_ring(ZZ)

    if n == 0:
        return Sequence([A], cr=True)

    F6_squared = F6**2
    D = delta_qexp(prec,var=var)
    Fprod = F6_squared
    Dprod = D
    R = A.parent()

    if cusp_only:
        ls = [R(0)] + [A] * n
        start = 1
    else:
        ls = [A] * (n+1)
        start = 0

    for i in range(1,n+1):
        ls[n-i] *= Fprod
        ls[i] *= Dprod

        Fprod *= F6_squared
        Dprod *= D
        Dprod = Dprod.add_bigoh(prec)

    if cusp_only:
        M = Matrix(ZZ, n, prec, [x.list() for x in ls[1:]])
        for i in range(n):
            for j in range(i):
                M.add_multiple_of_row(j,i,-M[j][i+1])
    else:
        M = Matrix(ZZ, n+1, prec, [x.list() for x in ls])
        for i in range(1,n+1):
            for j in range(i):
                M.add_multiple_of_row(j,i,-M[j][i])

    return Sequence([ R(x.list()).add_bigoh(prec) for x in M.rows() ], cr=True)

def delta_qexp(prec=10, var='q', K=ZZ):
    """
    Return the q-expansion of Delta as a power series with
    coefficients in K (=ZZ by default).

    INPUT:

        - ``prec`` -- integer; the absolute precision of the output

        - ``var`` -- (default: 'q') variable name

        - ``K`` -- (default: ZZ) base ring of answer

    OUTPUT:
        a power series over K


    ALGORITHM:
        Compute a simple very explicit modular form whose 8th power
        is Delta.   Then compute the 8th power using NTL polynomial
        arithmetic, which is VERY fast.   This function
        computes a *million* terms of Delta in under a minute.

    EXAMPLES::

        sage: delta_qexp(7)
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 + O(q^7)
        sage: delta_qexp(7,'z')
        z - 24*z^2 + 252*z^3 - 1472*z^4 + 4830*z^5 - 6048*z^6 + O(z^7)
        sage: delta_qexp(-3)
        Traceback (most recent call last):
        ...
        ValueError: prec must be positive

    AUTHORS:

        - William Stein: original code

        - David Harvey (2007-05): sped up first squaring step
    """
    if prec <= 0:
        raise ValueError, "prec must be positive"
    v = [0] * prec

    # Let F = \sum_{n >= 0} (-1)^n (2n+1) q^(floor(n(n+1)/2)).
    # Then delta is F^8.

    # First compute F^2 directly by naive polynomial multiplication,
    # since F is very sparse.

    stop = int((-1+math.sqrt(1+8*prec))/2.0)
    # make list of index/value pairs for the sparse poly
    values = [(n*(n+1)//2, ((-2*n-1) if (n & 1) else (2*n+1))) \
              for n in range(stop+1)]

    for (i1, v1) in values:
        for (i2, v2) in values:
            try:
                v[i1 + i2] += v1 * v2
            except IndexError:
                break

    # then use NTL to compute the remaining fourth power
    f = ntl.ZZX(v)
    t = verbose('made series')
    f = (f*f).truncate(prec)
    t = verbose('squared (2 of 3)', t)
    f = (f*f).truncate(prec)
    t = verbose('squared (3 of 3)', t)
    f = ntl.ZZX([0,1])*f
    t = verbose('shifted', t)
    R = K[[var]]
    f = R(f.truncate(prec), prec, check=False)
    t = verbose('coerced', t)
    return f

