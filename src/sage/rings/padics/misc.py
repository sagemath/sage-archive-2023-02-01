r"""
Miscellaneous Functions

This file contains two miscellaneous functions used by `p`-adics.

- ``min`` -- a version of ``min`` that returns `\infty` on empty input.
- ``max`` -- a version of ``max`` that returns `-\infty` on empty input.

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from six.moves.builtins import min as python_min
from six.moves.builtins import max as python_max
from sage.rings.infinity import infinity
from sage.rings.padics.factory import Zp
from sage.rings.all import PolynomialRing

def gauss_sum(a, p, f, prec = 20):
    """
    Returns the gauss sum g_q(a) = \sum_{u\in F_q^\times} omega(u)^(-a) \zeta_q^u
    where q = p^f, \omega is the Teichmuller character and \zeta_q is some arbitrary 
    choice of primitive p-th root of unity
    """
    a = a % (p**f)
    R = Zp(p, prec)
    R_poly = PolynomialRing(R,name='X')
    X = R_poly.gen()
    F = R.ext(X**(p-1)+p, names='pi')
    pi = F.gen()
    digits = Zp(p)(a).list(start_val = 0)
    n = len(digits)
    digits = digits + [0]*(f-n)
    s = sum(digits)
    out = -pi**(s)
    for i in range(0,f):
        a_i = R(sum([digits[k]*p**((i+k)%f) for k in range(f)]))        #an O(p^prec) term is necessay
        if a_i:
            out = out*R((a_i/(p**f-1)).gamma())                                #for coercing 0 correctly
    return out

def min(*L):
    r"""
    Return the minimum of the inputs, where the minimum of the empty
    list is `\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import min
        sage: min()
        +Infinity
        sage: min(2,3)
        2
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_min(L)
    except ValueError:
        return infinity


def max(*L):
    r"""
    Return the maximum of the inputs, where the maximum of the empty
    list is `-\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import max
        sage: max()
        -Infinity
        sage: max(2,3)
        3
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_max(L)
    except ValueError:
        return -infinity
