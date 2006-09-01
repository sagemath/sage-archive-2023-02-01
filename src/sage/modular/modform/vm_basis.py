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

from sage.matrix.all import MatrixSpace
from sage.modular.dims import dimension_cusp_forms_gamma0
from sage.rings.all import QQ, ZZ, Integer, binomial
from sage.structure.all import Sequence
from sage.libs.all import ntl
from sage.misc.all import verbose

from eis_series import eisenstein_series_qexp



def victor_miller_basis(k, prec=10, cusp_only=False, var='q'):
    r"""
    Compute and return the Victor-Miller basis for
    modular forms of weight k and level 1 to precision
    O(q^prec).  if \code{cusp_only} is True, return
    only a basis for the cuspidal subspace.

    INPUT:
        k -- an integer
        prec -- (default: 10) a positive integer
        cusp_only -- bool (default: False)
        var -- string (default: 'q'

    OUTPUT:
        list -- entries a power series in var defined over Q.

    EXAMPLES:
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
    """
    k = Integer(k)

    R = QQ[[var]]
    if k == 0:
        return Sequence([R(1, prec)], cr=True)
    elif k%2 == 1 or k < 4:
        return Sequence([])

    kk = k % 12
    if kk == 2:
        kk += 12
    b = None
    for a in range(15):
        c = kk - 4*a
        if c % 6 == 0:
            b = c // 6
            break
    assert not (b is None), "bug in VM basis"

    F4 = 240*eisenstein_series_qexp(4, prec=prec)
    F6 = -504*eisenstein_series_qexp(6, prec=prec)
    if var != 'q':
        F4 = R(F4)
        F6 = R(F6)
    Delta = (F4**3 - F6**2)/R(1728,prec)
    d = dimension_cusp_forms_gamma0(1, k)
    m = Delta / (F6*F6)
    g = m * F6**(2*d + b) * F4**a
    G = []
    for j in range(d):
        G.append(g)
        if j < d-1:
            g *= m

    if not cusp_only:
        G.insert(0, R(eisenstein_series_qexp(k, prec=prec)))

    M = MatrixSpace(QQ, len(G), prec)
    # we have to slice since precision in products can increase.
    e = [g.padded_list(prec) for g in G]
    A = M(sum(e, []))
    # this is still provably correct -- the guess is still proven right.
    # it's just that naive guess based on coefficients is way too big.
    E = A.echelon_form(height_guess=10**(k))
    return Sequence([R(list(v), prec) for v in E.rows()], cr=True)


## def delta_qexp(prec=10, var='q'):
##     """
##     Return the q-expansion of Delta.
##     """
##     F4 = 240*eisenstein_series_qexp(4, prec=prec)
##     F6 = -504*eisenstein_series_qexp(6, prec=prec)
##     R = QQ[[var]]
##     if var != 'q':
##         F4 = R(F4)
##         F6 = R(F6)
##     return (F4**3 - F6**2)/R(1728, prec)

def delta_qexp(prec=10, var='q'):
    """
    Return the q-expansion of Delta as a power series with
    coefficients in ZZ.

    ALGORITHM:
        Compute a simple very explicit modular form whose 8th power
        is Delta.   Then compute the 8th power using NTL polynomial
        arithmetic, which is VERY fast.   This function
        computes a *million* terms of Delta in under a minute.

    EXAMPLES:
        sage: delta_qexp(7)
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + O(q^7)
        sage: delta_qexp(7,'z')
        z - 24*z^2 + 252*z^3 - 1472*z^4 + 4830*z^5 - 6048*z^6 - 16744*z^7 + O(z^7)
        sage: delta_qexp(-3)
        Traceback (most recent call last):
        ...
        ValueError: prec must be positive
    """
    if prec <= 0:
        raise ValueError, "prec must be positive"
    v = [0] * prec
    stop = int((-1+math.sqrt(1+8*prec))/2.0)
    for n in range(stop+1):
        k = 2*n+1
        if n % 2 != 0:
            k = -k
        try:
            v[n*(n+1)//2] = k
        except IndexError:
            break

    f = ntl.ZZX(v)
    t = verbose('made series')
    f = (f*f).truncate(prec)
    t = verbose('squared (1 of 3)', t)
    f = (f*f).truncate(prec)
    t = verbose('squared (2 of 3)', t)
    f = (f*f).truncate(prec)
    t = verbose('squared (3 of 3)', t)
    f = ntl.ZZX([0,1])*f
    t = verbose('shifted', t)
    R = ZZ[[var]]
    f = R(f, prec, check=False)
    t = verbose('coerced', t)
    return f

