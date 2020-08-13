"""PolyBoRi's interface to libpolybori*

This file makes interfaces to PolyBoRi's runtime libraries available in Python via sage.


AUTHOR:
    The PolyBoRi Team, 2007-2012

            Examples:

            >>> from brial.frontend import *
            >>> r=declare_ring(["x0","x1","x2","y0","y1","y2"], globals())
            >>> x0>x1
            True
            >>> x0>x1*x2
            True
            >>> y0>y1
            True
            >>> y0>y1*y2
            True

            >>> r = r.clone(ordering=dlex)
            >>> r(x0) > r(x1)
            True
            >>> r(x0) > r(x1*x2)
            False

            >>> r = r.clone(ordering=dp_asc)
            >>> r(x0) > r(x1)
            False
            >>> r(x0) > r(x1*x2)
            False

            >>> r = r.clone(ordering=block_dlex, blocks=[3])
            >>> r(x0) > r(x1)
            True
            >>> r(x0) > r(x1*x2)
            False
            >>> r(x0) > r(y0*y1*y2)
            True

            >>> r = r.clone(ordering=block_dp_asc)
            >>> r(x0) > r(x1)
            False
            >>> r(x0) > r(y0)
            False
            >>> r(x0) > r(x1*x2)
            False

            >>> r = r.clone(ordering=block_dp_asc, blocks=[3])
            >>> r(x0) > r(y0)
            True

            >>> r(x0) > r(y0*y1)
            True

            >>> r = r.clone(names=["z17", "z7"])
            >>> [r.variable(idx) for idx in xrange(3)]
            [z17, z7, x2]
            >>> r = r.clone(names="abcde")
            >>> [r.variable(idx) for idx in xrange(6)]
            [a, b, c, d, e, y2]

"""

from sage import all
from sage.rings.polynomial.pbori.pbori.pbori import *

import weakref


OrderCode = type('OrderCode', (object,), order_dict)


def Ring(n, order='lp', names=None, blocks=[]):

    pbnames = names
    if pbnames is None:
        pbnames = ['x(' + str(idx) + ')' for idx in range(n)]
    order = TermOrder_from_pb_order(n, order, blocks)
    R = BooleanPolynomialRing(n, names=pbnames, order=order)
    return R

BoolePolynomialVector = BooleanPolynomialVector


#todo: PolyBoRi's original interface uses its WeakRingPtr here
def WeakRingRef(ring):
    return weakref.weakref(ring)

Monomial = MonomialFactory()
Polynomial = PolynomialFactory()
Variable = VariableFactory()


_add_up_polynomials = add_up_polynomials


def add_up_polynomials(polys, init):
    """
    Adds up the polynomials in polys (which should be a BoolePolynomialVector or a sequence of ???
    """
    if not isinstance(polys, BoolePolynomialVector):
        vec = BoolePolynomialVector
        for p in polys:
            vec.append(p)
        polys = vec

    return _add_up_polynomials(polys, init)

_gauss_on_polys = gauss_on_polys


def gauss_on_polys(l):
    vec = BoolePolynomialVector(l)
    return list(_gauss_on_polys(vec))
